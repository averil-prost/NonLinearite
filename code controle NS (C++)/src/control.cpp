/**
 * Résolution de Navier-Stokes instationnaire + contrôle
 * PFE/mémoire MFA 2021-2022, département GM INSA Rouen Normandie
 * Averil PROST (averil.prost@insa-rouen.fr)
 * sous la direction d'Antoine Tonnoir & Nicolas Forcadel.
 */

#include "../include/datastructure.hpp"
#include "../include/visu.hpp"
#include "../include/procedure/assemblage.hpp"
#include "../include/procedure/evaluation.hpp"
#include "../include/procedure/iterations.hpp"
#include "../include/interface.hpp"

#define CFL_CONSTANT 5.0
#define VISCOSITY 0.1 // 5e-5 
#define SCHEMA "navierim"

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
////
////    MAIN
////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

int main () {
    std::cout << "Welcome in insta." << std::endl;

    ////////////////////////////////////////////////////////////////
    // Paramètres
    ////////////////////////////////////////////////////////////////

    // contrôle d'exécution
    std::stringstream runname_root; runname_root << "gunz"; // nom qui sera associé aux fichiers de sortie .vtk
    bool load_matrices       = true;        // l'existence des fichiers n'est pas vérifiée
    bool save_matrices       = true;        // pas d'effet si load_matrices est à true
    bool export_mesh         = !true;       // uniquement pour les maillages de la validation
    bool iterate             = true;        
    bool validate            = !true;       
    bool write_vtk           = true;        
    bool verbose             = true;        
    bool checkStationnary    = !true;       // si la solution PEUT se stabiliser (évite NaN, mais + long)

    // maillage
    // std::string mesh_name_root = "square";
    std::string mesh_name_root = "fly_coarse";
    std::vector<double> mesh_scales = {0.15};//{0.25, 0.15, 0.1};

    // schéma
    std::string schema_name = SCHEMA; // stokes, navierex, navierim
    int deg_ef_v = 2;
    int deg_ef_p = 1;

    // simulation
    double nu = VISCOSITY;
    double T = 1.0; 
    int Nit_max = 20; // nombre de pas de temps calculés
    int Npass_max = 100; // nombre de cycles forward-backward maximum
    int Ncutstep_max = 10; // nombre de réductions du coût max
    double epsilon = 1e-5; // stopping criteria on J norm

    // paramètres de la fonction coût
    double alpha = 1.0; // pondère la condition sur tout [0,T]
    double beta = 0.0001; // pondère la norme du contrôle
    double gamma = 0.5; // pondère la condition à t=T

    // post-traitement
    std::vector<std::string> reports;
    std::vector<std::vector<double>> Jvalues, Fnorms, internErrs; // oui, j'aime la mémoire
    std::string deco = "///////////////////////////////////////////////////////////////";

    runname_root << (schema_name=="stokes"?"S":(schema_name=="navierex"?"NE":"NI")) \
            << "_b"<< (int)1e4*beta << "_g" << (int)1e2*gamma << "_v" << (int)1e2*VISCOSITY;

    for (double mesh_scale : mesh_scales) {
        std::stringstream report;
        std::vector<double> Jvalue, Fnorm, internErr, finalErr;

        ////////////////////////////////////////////////////////////////
        // Maillage
        ////////////////////////////////////////////////////////////////

        // pour ce code, mesh_scale > 0.01 => on évite les . dans les noms de fichiers
        std::string mesh_name = mesh_name_root+(export_mesh?std::to_string((int)(100*mesh_scale)):"");
        if (export_mesh) {
            interface::exportSquareMesh (MESH_FOLDER+mesh_name, mesh_scale);
        }

        // récupérer le maillage (indispensable pour pouvoir écrire les .vtk...)
        Meshes mesh, mesh_v, mesh_p;
        char name [40]; sprintf(name, "%s", (MESH_FOLDER+mesh_name + ".msh").c_str());
        mesh.readGMSH (name); // je sais que ce n'est pas propre de le lire 3 fois
        mesh_v.readGMSH (name); // mais après beaucoup de galère avec la mémoire, 
        mesh_p.readGMSH (name); // c'est finalement bien plus simple.
        double h = mesh.getEleSize2D(); // estimation de la mesure d'un triangle

        element_fini_D2 ef_v (deg_ef_v), ef_p (deg_ef_p);
        mesh_v.convertP1ToPk2D (ef_v.getTaille(), ef_v.getTreilli()); // modifie irréversiblement !!!!
        mesh_p.convertP1ToPk2D (ef_p.getTaille(), ef_p.getTreilli()); // modifie irréversiblement !!!!

        ////////////////////////////////////////////////////////////////
        // Données diverses du schéma
        ////////////////////////////////////////////////////////////////

        // construction d'un ~ gestionnaire de bordure ~ (obligatoire... mais ça va, pas très long)
        datastructure::bordermgr bman = builder::bordermgr (mesh_v, ef_v.getTaille());

        // Dirichlet géré par pseudo-élimination
        int cardBaseV = mesh_v.getNNodes(); 
        int cardBaseP = mesh_p.getNNodes();

        // condition CFL : stokes et narvierim inconditionnellement stables
        int Nit=Nit_max;
        double deltaT = 1.0*T/Nit;
        if (schema_name == "navierex") {
            deltaT = CFL_CONSTANT * nu * h*h;
            Nit = std::ceil (T/deltaT);
        }

        // assemblage des matrices des schémas
        datastructure::schemas schemas = procedure::assemblerschemas (schema_name, mesh_name, ef_v, ef_p, 
            load_matrices, save_matrices, mesh_v, mesh_p, bman, nu, deltaT, cardBaseV, cardBaseP);

        report << deco << "\nProgramme de " << Nit << " itérations, cardBaseV/P : " << cardBaseV << "/" << cardBaseP << std::endl;
        if (verbose) std::cout << deco << "\nProgramme de " << Nit << " itérations, cardBaseV/P : " << cardBaseV << "/" << cardBaseP << std::endl;

        if (iterate) {
            ////////////////////////////////////////////////////////////////
            // Main initialisation
            ////////////////////////////////////////////////////////////////

            int ipass=0, icutstep=0, n=0;
            bool stopPass=false, stopCutStep=false;
            double costn=0.0, costnp1=0.0, gradstep=1.0, internError=0.0, finalError=0.0, controlNorm=0.0;

            // all variables goes from 0 to Nit
            // objectif et relèvements dudit objectif
            std::vector<Vecteur> U_Q, U_Qborder;
            // respectively state, adjoint state and control. 
            std::vector<Vecteur> U, W, F, prevF; 

            // initialisation nulle pour F, prevF et W
            for (n=0; n<=Nit; ++n) {
                Vecteur theW (2*cardBaseV+cardBaseP+1), theF (2*cardBaseV+cardBaseP+1), theprevF (2*cardBaseV+cardBaseP+1);
                W.push_back (theW); F.push_back (theF); prevF.push_back (theprevF);
            }

            // interpoler la cible et avoir ses relèvements
            procedure::buildTarget (U_Q, U_Qborder, procedure::U_Q, mesh_v, Nit, deltaT, cardBaseV, cardBaseP, bman);

            // cible : solution du problème de Stokes stationnaire
            // procedure::buildTarget (U_Q, U_Qborder, procedure::dirichlet, nu, mesh_v, schemas, Nit, cardBaseV, cardBaseP, bman);

            // // résolution du problème stationnaire pour la condition initiale (faut bien avoir quelque chose...)
            // Vecteur initialScdMembre = procedure::initialScdMembre (schemas.insta.uv, schemas.insta.gradugradv,  \
            //     schemas.insta.mdivu1v2, nu, F.at(0), U_Qborder.at(0), bman, cardBaseV, cardBaseP);
            // U.push_back(schemas.sta.tobeinverted.solveIte (initialScdMembre, U_Q.at(0), "KrylovSym", "Id", 1.47, 8000, 1000, 1e-8, 0));

            // initialisation de U l'état
            for (n=0; n<=Nit; ++n) {
                U.push_back (U_Q.at(n) * (-10.0)); // following example 1, Gunzburger : u0 = -10 * U_Q
                // U.push_back (U_Q.at(n) * (0.0)); 
            }

            std::string entetev_vtk = interface::entete_vtk (mesh_v.getNodes(), mesh_v.convertPkToP12D(ef_v.getTaille(), ef_v.getTreilli(), ef_v.getDecoup()));
            std::string entetep_vtk = interface::entete_vtk (mesh_p.getNodes(), mesh_p.convertPkToP12D(ef_p.getTaille(), ef_p.getTreilli(), ef_p.getDecoup()));
            interface::writeVTKAtStepN ("sol_"+mesh_name_root, n, entetev_vtk, entetep_vtk, U_Q.at(0), U_Qborder.at(0)*0.0, cardBaseV, cardBaseP);

            // first forward pass
            internError = procedure::norme_v_carre (U.at(0), U_Q.at(0), cardBaseV); 
            controlNorm = procedure::norme_v_carre (F.at(0), cardBaseV);
            if (verbose) std::cout << "//// Initial forward pass" << std::endl;
            for (n=0; n<Nit; ++n) {
                U.at(n+1) = procedure::iterate_forward (schemas.insta, U.at(n), F.at(n+1), U_Qborder.at(n+1), U_Qborder.at(n), bman, deltaT, nu, checkStationnary, verbose);

                internError += procedure::norme_v_carre (U.at(n+1), U_Q.at(n+1), cardBaseV); 
                controlNorm += procedure::norme_v_carre (F.at(n+1), cardBaseV);
            } // for n forward
            // fin de calcul du coût
            finalError = procedure::norme_v_carre (U.at(Nit), U_Q.at(Nit), cardBaseV);
            costn = deltaT * (alpha * internError + beta * controlNorm) + gamma * 0.5 * finalError;

            ////////////////////////////////////////////////////////////////
            // Main loop
            ////////////////////////////////////////////////////////////////

            while (ipass < Npass_max && !stopPass) {
                if (verbose) std::cout << deco+"\n" << "//   Pass " << ipass << "\n"+deco << std::endl;

                /////////// Backward pass
                W.at(Nit) = (U.at(Nit) - U_Q.at(Nit)) * gamma;
                // backward iterations
                if (verbose) std::cout << "//// Backward pass " << ipass << std::endl;
                for (n=Nit-1; n>=0; --n) {
                    W.at(n) = procedure::iterate_backward (schemas.insta, W.at(n+1), U.at(n+1), U_Q.at(n+1), bman, deltaT, verbose);

                    // std::cout << visu::tostring (W.at(n), cardBaseV, cardBaseP);
                } // for n backward

                stopCutStep = false; icutstep=0;
                while (!stopCutStep && icutstep < Ncutstep_max) { // reducing gradstep until we diminish the cost
                    /////////// Control update 
                    for (n=0; n<=Nit; ++n) {
                        F.at(n) = prevF.at(n) * (1.0 - gradstep * beta) - W.at(n) * gradstep;
                    }

                    /////////// Forward pass + cost evaluation
                    // u0 fixé au départ

                    internError = procedure::norme_v_carre (U.at(0), U_Q.at(0), cardBaseV); 
                    controlNorm = procedure::norme_v_carre (F.at(0), cardBaseV);
                    if (verbose) std::cout << "//// Forward pass " << ipass << std::endl;
                    for (n=0; n<Nit; ++n) {
                        U.at(n+1) = procedure::iterate_forward (schemas.insta, U.at(n), F.at(n), U_Qborder.at(n+1), U_Qborder.at(n), bman, deltaT, nu, checkStationnary, verbose);

                        internError += procedure::norme_v_carre (U.at(n+1), U_Q.at(n+1), cardBaseV); 
                        controlNorm += procedure::norme_v_carre (F.at(n+1), cardBaseV);
                    } // for n forward
                    // fin de calcul du coût
                    finalError = procedure::norme_v_carre (U.at(Nit), U_Q.at(Nit), cardBaseV);
                    costnp1 = deltaT * (alpha * internError + beta * controlNorm) + gamma * 0.5 * finalError;

                    /////////// stopping criteria
                    if (costnp1 > costn) {
                        gradstep *= 0.5;
                        ++icutstep;
                        std::cout << "\t\t\t\t\tCUTTING" << std::endl;
                    }
                    else if (std::fabs(costnp1 - costn) / costn > epsilon) {
                        gradstep *= 1.5;
                        stopCutStep = true;
                    }
                    else {
                        stopPass = true;
                        stopCutStep = true;
                    }
                } // while cutStep

                costn = costnp1;
                prevF = F;
                ++ipass;

                report << "Valeur de J à l'étape " << ipass << " : " << std::setprecision(10) << costn \
                    << ", intern : " << internError << ", final : " << finalError << ", Fnorm : " << controlNorm << std::endl;
                if (verbose) std::cout << "/////\t\tValeur de J à l'étape " << ipass << " : " << std::setprecision(10) << costn \
                    << ", intern : " << internError << ", final : " << finalError << ", Fnorm : " << controlNorm << std::endl;
                Jvalue.push_back (costn);
            } // while pass

            ////////////////////////////////////////////////////////////////
            // Enregistrement des résultats
            ////////////////////////////////////////////////////////////////

            std::string runname = runname_root.str() + mesh_name;
            if (write_vtk) {
                std::string entetev_vtk = interface::entete_vtk (mesh_v.getNodes(), mesh_v.convertPkToP12D(ef_v.getTaille(), ef_v.getTreilli(), ef_v.getDecoup()));
                std::string entetep_vtk = interface::entete_vtk (mesh_p.getNodes(), mesh_p.convertPkToP12D(ef_p.getTaille(), ef_p.getTreilli(), ef_p.getDecoup()));
                for (n=0; n<Nit; ++n) {
                    interface::writeVTKAtStepN (runname, n, entetev_vtk, entetep_vtk, U.at(n), U_Qborder.at(n), cardBaseV, cardBaseP);
                    interface::writeVTKAtStepN (runname+"_control", n, entetev_vtk, entetep_vtk, F.at(n), U_Qborder.at(n)*0.0, cardBaseV, cardBaseP);
                }
            } // writing results

            for (n=0; n<=Nit; ++n) {
                Fnorm.push_back (std::sqrt(procedure::norme_v_carre(F.at(n), cardBaseV))); 
                internErr.push_back (std::sqrt(procedure::norme_v_carre(U.at(n), U_Q.at(n), cardBaseV))); 
            }

            Jvalues.push_back (Jvalue);
            Fnorms.push_back (Fnorm);
            internErrs.push_back (internErr);

        } // if iterate

        reports.push_back(report.str());
    } // for meshes


    ////////////////////////////////////////////////////////////////
    // Post-traitement
    ////////////////////////////////////////////////////////////////

    // affichage des rapports 
    for (std::string report : reports) {
        std::cout << report;
    }

    // enregistrement des mesures d'exécution
    interface::saveMeasure<double> (Jvalues, ERR_FOLDER+runname_root.str()+"_Jvalue.txt");
    interface::saveMeasure<double> (Fnorms,  ERR_FOLDER+runname_root.str()+"_Fnorm.txt");
    interface::saveMeasure<double> (internErrs,  ERR_FOLDER+runname_root.str()+"_internErr.txt");

    std::cout << "Bye" << std::endl;
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////