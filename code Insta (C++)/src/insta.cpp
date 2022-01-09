/**
 * Résolution de Navier-Stokes instationnaire
 * PFE 2021-2022, département GM INSA Rouen Normandie
 * Averil PROST (averil.prost@insa-rouen.fr)
 * sous la direction d'Antoine Tonnoir.
 */

#include "../include/datastructure.hpp"
#include "../include/visu.hpp"
#include "../include/procedure/assemblage.hpp"
#include "../include/procedure/iterations.hpp"
#include "../include/procedure/test.hpp"
#include "../include/interface.hpp"

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
////
////    DATA
////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

#define CFL_CONSTANT 5.0
#define VISCOSITY 0.1
#define EXAMPLE 3
#define SCHEMA "navierim"

// terme source
Vecteur f (Vecteur& x, double t) {
    Vecteur res = x*0.0; // default

    #if (EXAMPLE==2)
        if (SCHEMA != "stokes") {
            res.set(0, - 6.0*x(1)); // seulement si navierim/ex
        }
        res.set(1, 1.0);
    #endif

    #if (EXAMPLE==3)
        if (t>1.0 && t<=2.0) {
            res.set(0, - x(1)*(1.0-x(1))*M_PI*0.5*std::sin(M_PI*0.5*(t-1.0))); 
            res.set(1,   x(0)*(1.0-x(0))*M_PI*0.5*std::cos(M_PI*0.5*(t-1.0))); 
            if (SCHEMA!="stokes") {
                res.set(0, res(0) + std::sin(M_PI*0.5*(t-1.0))*std::cos(M_PI*0.5*(t-1.0))*(1.0-2.0*x(1))*x(0)*(1.0-x(0))); 
                res.set(1, res(1) + std::sin(M_PI*0.5*(t-1.0))*std::cos(M_PI*0.5*(t-1.0))*(1.0-2.0*x(0))*x(1)*(1.0-x(1)));
            }
        }
    #endif 

    return res;
}

// coord x de la vitesse
double solv1 (Vecteur& x, double t) {
    #if (EXAMPLE==1)    
        return x(1)*(1.0-x(1)); // square, example 1
    #endif 
    #if (EXAMPLE==2)    
        return -3.0*x(1)*x(1); // square, example 2
    #endif 
    #if (EXAMPLE==3)    
        return x(1)*(1.0-x(1)) * (t<=1.0 ? 1.0 : (t <= 2.0 ? std::cos(M_PI*0.5*(t-1.0)) : 0.0));
    #endif 
}
// coord y de la vitesse
double solv2 (Vecteur& x, double t) {
    #if (EXAMPLE==1)    
        return 0.0; // square, example 1
    #endif 
    #if (EXAMPLE==2)    
        return 1.0; // square, example 2
    #endif 
    #if (EXAMPLE==3)    
        return x(0)*(1.0-x(0)) * (t<=1.0 ? 0.0 : (t <= 2.0 ? std::sin(M_PI*0.5*(t-1.0)) : 1.0));
    #endif 
}
// pression
double solp (Vecteur& x, double t) {
    #if (EXAMPLE==1)    
        return - 2.0 * VISCOSITY * x(0) + VISCOSITY; // square, example 1
    #endif 
    #if (EXAMPLE==2)    
        return x(1) - 6.0 * VISCOSITY * x(0) - (0.5 - 3.0*VISCOSITY); // square, example 2
    #endif 
    #if (EXAMPLE==3)    
        double res = 0.0;
        if (t <= 1.0) {
            res = - 2.0 * VISCOSITY * x(0) + VISCOSITY;
        } else if (t <= 2.0) {
            res = (-2.0 * VISCOSITY * x(0) + VISCOSITY)*std::cos(M_PI*0.5*(t-1.0)) \
                + (-2.0 * VISCOSITY * x(1) + VISCOSITY)*std::sin(M_PI*0.5*(t-1.0));
        }
        else {
            res = - 2.0 * VISCOSITY * x(1) + VISCOSITY;
        }
        return res;
    #endif 
}

Vecteur dirichlet (Vecteur& x, double t) {
    Vecteur res (2);
    // pour la validation
    #if (EXAMPLE <= 3)
    res.set(0, solv1(x,t));
    res.set(1, solv2(x,t));
    #endif 

    #if (EXAMPLE == 4) // pour l'exemple arrow
    if (!(x(0) >= 0.4 && x(0) <= 1.0 && x(1) >= 0.1 && x(1) <= 2.4)) {
        // autour de l'obstacle : nada. Sinon, profil de Poiseuille
        res.set(0, x(1)*(2.0-x(1)));
    }
    #endif 

    #if (EXAMPLE == 5) // pour karman
    if (!(x(0) >= 0.1 && x(0) <= 1.0 && x(1) >= 0.2 && x(1) <= 0.8)) {
        // autour de l'obstacle : nada. Sinon, profil de Poiseuille
        res.set(0, 0.05);
        if (x(0) < 0.1) { // pour les proches de 0 : petite perturbation (inutile)
            res.set(1, 0.05*(0.1 - x(0))*std::cos(0.5*M_PI*t)*x(1)*(1.0-x(1)));
        }
    } 
    #endif
    return res;
}

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
    std::string runname_root = "ex3";   // nom qui sera associé aux fichiers de sortie .vtk
    bool load_matrices       = !true;        // l'existence des fichiers n'est pas vérifiée
    bool save_matrices       = true;        // pas d'effet si load_matrices est à true
    bool export_mesh         = true;       // uniquement pour les maillages de la validation
    bool iterate             = true;        
    bool validate            = true;       
    bool write_vtk           = true;        
    bool verbose             = true;        
    bool stationary_f        = !true;        // terme source 
    bool stationary_d        = !true;       // dirichlet
    bool checkStationnary    = true;       // si la solution PEUT se stabiliser (évite NaN, mais + long)

    // maillage
    std::string mesh_name_root = "square";
    std::vector<double> mesh_scales = {0.2, 0.15, 0.1, 0.09, 0.08}; // pour la validation

    // schéma
    std::string schema_name = SCHEMA; // stokes, navierex, navierim
    int deg_ef_v = 2;
    int deg_ef_p = 1;

    // physique (enfin le plus proche)
    double nu = VISCOSITY; // viscosité, garder dans les [e-5,e-1]
    double T = 3.0; // 10 frame par seconde
    int Nit = 30;
    // valeurs pour Karman : viscosité 5e-5, T = 200.0, Nit = 2000

    // post-traitement
    std::vector<std::string> reports, errors;

    for (double mesh_scale : mesh_scales) {

        ////////////////////////////////////////////////////////////////
        // Maillage
        ////////////////////////////////////////////////////////////////

        // pour ce code, mesh_scale > 0.01 => on évite les . dans les noms de fichiers (validation)
        std::string mesh_name = mesh_name_root;
        if (validate) mesh_name += std::to_string((int)(100*mesh_scale));
        if (validate && export_mesh) {
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
        if (verbose) std::cout << "cardBaseV : " << cardBaseV << ", cardBaseP : " << cardBaseP << std::endl;

        // condition CFL : stokes et narvierim inconditionnellement stables
        double deltaT = 1.0;
        if (schema_name == "navierex") {
            deltaT = CFL_CONSTANT * nu * h*h; // ! au péril de l'utilisateur
            Nit = std::ceil (T/deltaT);
        }
        deltaT = T/Nit;
        if (verbose) std::cout << "Programme de " << Nit << " itérations avec deltaT=" << deltaT << " et T=" << T << std::endl;

        // préparation au calcul d'une solution de référence (pour l'erreur de quadrature)
        std::vector<int> correspondance_v, correspondance_p;
        if (validate) {
            correspondance_v = builder::get_corresp (ef_v, mesh_v);
            correspondance_p = builder::get_corresp (ef_p, mesh_p);
        }

        // note au lecteur : ce n'est malheureusement pas parce que c'est caché dans une fonction que c'est plus joli
        datastructure::schemas schemas = procedure::assemblerschemas (schema_name, mesh_name, ef_v, ef_p, 
            load_matrices, save_matrices, mesh_v, mesh_p, bman, nu, deltaT, cardBaseV, cardBaseP);
        if (verbose) std::cout << "Schéma construit." << std::endl;

        if (iterate) {
            ////////////////////////////////////////////////////////////////
            // Initialisation
            ////////////////////////////////////////////////////////////////

            Vecteur fn (2*cardBaseV+cardBaseP+1), dirich1 (2*cardBaseV+cardBaseP+1), dirich2 (2*cardBaseV+cardBaseP+1);
            // relèvement : on a juste besoin de celui à n et n-1, donc alterne dans le stockage (courant = n%2)
            std::vector<Vecteur> dirich; dirich.push_back(dirich1); dirich.push_back(dirich2);
            Vecteur Un (2*cardBaseV+cardBaseP+1), Unp1 (2*cardBaseV+cardBaseP+1);

            // interpoler le terme source 
            builder::interpolate (fn, f, mesh_v, 0.0);  
            // construire la donnée initiale
            if (validate) {
                #if (EXAMPLE != 1 && EXAMPLE != 2) // sinon on part de u_0 = 0
                    Un = builder::getSol_interp (solv1, solv2, solp, 0.0, mesh_v, mesh_p, cardBaseV, cardBaseP);
                #endif 

                if (verbose) std::cout << visu::stats (Un, cardBaseV, cardBaseP, 
                    procedure::erreurL2_v(Un,dirich.at(0),solv1,solv2,ef_v,correspondance_v,mesh_v,cardBaseV,0.0,bman),
                    procedure::erreurL2_p (Un,solp,ef_p,correspondance_p,mesh_p,cardBaseV,0.0), 
                    procedure::test::checkPression (Un,schemas.insta.tobeinverted,cardBaseV,cardBaseP));
            } else {
                // résolution du problème stokes stationnaire pour la condition initiale
                Vecteur initialScdMembre = procedure::initialScdMembre (schemas.insta.uv, schemas.insta.gradugradv,  \
                    schemas.insta.mdivu1v2, nu, fn, dirich.at(0), bman, cardBaseV, cardBaseP);
                Un = schemas.sta.tobeinverted.solveIte (initialScdMembre, fn, "KrylovSym", "Id", 1.47, 8000, 1000, 1e-8, 0);
            }

            // gérer le relèvement
            if (stationary_d) {
                builder::relever(dirich, dirichlet, mesh_v, 0, 0.0, cardBaseV, bman);
                dirich.at(1) = dirich.at(0);
            }
            // projection de la condition initiale : on élimine les parasites au bord
            for (int inode = 0; inode < cardBaseV; ++inode) {
                if (!bman.isInside.at(inode)) {
                    if (!stationary_d) { // uniquement pour la première écriture, sera écrasé
                        dirich.at(0).set(2*inode  , Un(2*inode  )); 
                        dirich.at(0).set(2*inode+1, Un(2*inode+1)); 
                    }
                    Un.set(2*inode  , 0.0); 
                    Un.set(2*inode+1, 0.0); 
                }
            }

            // préparation à l'écriture des fichiers de sortie
            std::stringstream suffix; suffix << "ev" << deg_ef_v << "_ep" << deg_ef_p << "nu" << (int)1e5*VISCOSITY;
            std::string runname = runname_root+(schema_name=="stokes"?"S":(schema_name=="navierex"?"NE":"NI")) + mesh_name + suffix.str();
            std::string entetev_vtk, entetep_vtk; // à faire après que les mesh_v et mesh_p aient été convertis
            if (write_vtk) {
                entetev_vtk = interface::entete_vtk (mesh_v.getNodes(), mesh_v.convertPkToP12D(ef_v.getTaille(), ef_v.getTreilli(), ef_v.getDecoup()));
                entetep_vtk = interface::entete_vtk (mesh_p.getNodes(), mesh_p.convertPkToP12D(ef_p.getTaille(), ef_p.getTreilli(), ef_p.getDecoup()));
                // écriture de l'état initial 
                interface::writeVTKAtStepN (runname, 0, entetev_vtk, entetep_vtk, Un, dirich.at(0), cardBaseV, cardBaseP); 
            }      

            ////////////////////////////////////////////////////////////////
            // Itérations 
            ////////////////////////////////////////////////////////////////

            double erreur_v = 0.0, erreur_p=0.0, intP=0.0, erreur_v_totale=0.0, erreur_p_totale=0.0; // validation
            for (int n=1; n<=Nit; ++n) {
                if (verbose) std::cout << "/////////////////////////// Itération " << n << " / " << Nit << " ///////////////////////////" << std::endl;

                // itérations elles-mêmes
                if (!stationary_d) builder::relever(dirich, dirichlet, mesh_v, n%2, n*deltaT, cardBaseV, bman);
                if (!stationary_f) builder::interpolate (fn, f, mesh_v, n*deltaT); 
                Unp1 = procedure::iterate (schemas.insta, Un, fn, dirich.at(n%2), dirich.at((n-1)%2), bman, n, deltaT, nu, checkStationnary); 
                Un = Unp1;

                // trucs fun en plus : validation, écriture des erreurs et des résultats
                if (validate) {
                    erreur_v = procedure::erreurL2_v (Unp1, dirich.at(n%2), solv1, solv2, ef_v, correspondance_v, mesh_v, cardBaseV, n*deltaT, bman);
                    erreur_p = procedure::erreurL2_p (Unp1, solp, ef_p, correspondance_p, mesh_p, cardBaseV, n*deltaT);
                    erreur_v_totale += erreur_v; erreur_p_totale += erreur_p;
                    intP = procedure::test::checkPression (Unp1, schemas.insta.tobeinverted, cardBaseV, cardBaseP);
                    if (verbose) std::cout << visu::stats (Unp1, cardBaseV, cardBaseP, erreur_v, erreur_p, intP);
                }
                if (write_vtk) interface::writeVTKAtStepN (runname, n, entetev_vtk, entetep_vtk, Un, dirich.at(n%2), cardBaseV, cardBaseP);  
            } // for n

            std::cout << "//////////////////////////////////////////////////////" << std::endl;
            if (validate) {
                // pour tout afficher à la fin, après les verbose
                std::stringstream report;
                report << "Maillage " << mesh_name << " à " << mesh.getNNodes() << " noeuds\t: ";
                report << std::fixed << std::setprecision(16);
                report << "erreur vitesse " << erreur_v_totale*deltaT << ", pression " << erreur_p_totale*deltaT << std::endl;
                reports.push_back (report.str());
                // pour les fichiers de convergence
                std::stringstream error;
                error << std::fixed << std::setprecision(16);
                error << mesh.getNNodes() << "\t";
                error << h << "\t";
                error << erreur_v_totale*deltaT << "\t";
                error << erreur_p_totale*deltaT << std::endl;
                errors.push_back(error.str());
            } // if validate
        } // if iterate
    } // for meshes

    // affichage des rapports 
    std::cout << runname_root << std::endl;
    for (auto report : reports) {
        std::cout << report;
    }

    if (validate) {
        // enregistrement des erreurs
        std::ofstream errstream ("data/erreurs/square_"+schema_name+"_ex"+std::to_string(EXAMPLE)+".txt");
        errstream << "#nombre de points - h - erreur vitesse - erreur pression" << std::endl;
        for (auto error : errors) {
            errstream << error;
        }
        errstream.close ();
    }    

    std::cout << "Bye" << std::endl;
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////