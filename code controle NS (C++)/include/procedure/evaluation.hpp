#ifndef __EVALUATION__
#define __EVALUATION__

#include "../inclusions.hpp"
#include "../datastructure.hpp"
#include "../builder.hpp"
#include "../interface.hpp"
#include "../visu.hpp"

// following example 1, p. 1503 (24 on pdf) of Velocity tracking for Navier-Stokes flows
// Analysis and approximation of [...] with distributed controls, Gunzburger & Manservisi
double phi (double z) {
    return (1.0 - std::cos(0.8*M_PI*z))*(1.0-z)*(1.0-z);
}

double phiPrime (double z) {
    return (1.0 - z)*(std::sin(0.8*M_PI*z)*(1-z)*0.8*M_PI - 2.0*(1.0 - std::cos(0.8*M_PI*z)));
}

namespace procedure {

    // dirichlet imposé 
    Vecteur dirichlet (Vecteur& x) {
        Vecteur res (2);
        if (!((x(0) > 0.0) && (x(0) < 1.7) && (x(1) > 0.1) && (x(1) < 2.0))) {
            res.set(0, 1.0);
            res.set(1, -1.0);
        }
        return res;
    }

    // objectif à atteindre (en vitesse)
    Vecteur U_Q (Vecteur& x, double t) {
        Vecteur res (2);
        // gunzburger example
        res.set(0,   10.0 * phi(x(0)) * phiPrime (x(1)));
        res.set(1, - 10.0 * phi(x(1)) * phiPrime (x(0)));

        // Poiseuille
        // res.set (0, x(1)*(1.0-x(1)));
        return res;
    }


    // renvoie tous les U_Q ainsi que les relèvements
    void buildTarget (std::vector<Vecteur>& U_Q, std::vector<Vecteur>& U_Qborder, Vecteur source (Vecteur&, double),
        Meshes& mesh_v, int Nit, double deltaT, int cardBaseV, int cardBaseP, datastructure::bordermgr& bman) {
        std::cout << "Entrée dans buildTarget" << std::endl;
        int inode=0, n=0; 
        std::vector<Vecteur*> nodes = mesh_v.getNodes();
        Vecteur value (3);
        for (n=0; n<=Nit; ++n) { // pour chaque temps
            Vecteur target (2*cardBaseV+cardBaseP+1), targetborder (2*cardBaseV+cardBaseP+1);
            inode=0;
            for (auto node : nodes) { // pour chaque noeud du maillage
                value = source (*node, n*deltaT);
                target.set(2*inode  , value(0)); // interpolation
                target.set(2*inode+1, value(1));
                if (!bman.isInside.at(inode)) {
                    targetborder.set(2*inode  , value(0));
                    targetborder.set(2*inode+1, value(1));
                }
                ++inode;
            }
            U_Q.push_back (target);
            U_Qborder.push_back (targetborder);
        } // for n
        std::cout << "Sortie de buildTarget" << std::endl;
    }

    // construction de U_Q et U_Qborder pour une solution stationnaire du problème de Stokes
    void buildTarget (std::vector<Vecteur>& U_Q, std::vector<Vecteur>& U_Qborder, Vecteur dirich (Vecteur&), double nu,
        Meshes& mesh_v, datastructure::schemas& schemas, int Nit, int cardBaseV, int cardBaseP, datastructure::bordermgr& bman) {
        std::cout << "Entrée dans buildTarget par Stokes sta" << std::endl;

        // interpolation de la condition au bord
        Vecteur target (2*cardBaseV+cardBaseP+1), targetborder (2*cardBaseV+cardBaseP+1);
        Vecteur value (2); int inode=0;         
        for (auto node : mesh_v.getNodes()) { // pour chaque noeud du bord maillage
            if (!bman.isInside.at(inode)) {
                value = dirich (*node);
                targetborder.set(2*inode  , value(0));
                targetborder.set(2*inode+1, value(1));
            }
            ++inode;
        }
        // construction du second membre du système de Stokes
        // attention malversation : target n'a rien à faire là, c'est juste un vecteur nul
        Vecteur scdmembre = procedure::initialScdMembre (schemas.insta.uv, schemas.insta.gradugradv,  \
            schemas.insta.mdivu1v2, nu, target, targetborder, bman, cardBaseV, cardBaseP);
        // résolution une bonne fois pour toute
        target = schemas.sta.tobeinverted.solveIte (scdmembre, target, "GMRES", "Id", 1.47, 8000, 1000, 1e-8, 0);

        // forcer l'égalité au bord
        for (inode=0; inode<cardBaseV; ++inode) {
            if (!bman.isInside.at(inode)) {
                target.set(2*inode  , targetborder(2*inode  ));
                target.set(2*inode+1, targetborder(2*inode+1));
            }
        }

        // solution stationnaire 
        for (int n=0; n<=Nit; ++n) { // pour chaque temps
            Vecteur datarget = target, datargetborder = targetborder;
            U_Q.push_back (datarget);
            U_Qborder.push_back (datargetborder);
        } // for n
        std::cout << "Sortie de buildTarget par Stokes sta" << std::endl;
    }

    // carré de la norme de la partie vitesse
    double norme_v_carre (Vecteur& U, int cardBaseV) {
        double res;
        for (int i=0; i<cardBaseV; ++i) {
            res += U(2*i)*U(2*i) + U(2*i+1)*U(2*i+1);
        }
        return res;
    }

    // carré de la norme de la partie vitesse de U - V
    double norme_v_carre (Vecteur& U, Vecteur& V, int cardBaseV) {
        double res;
        for (int i=0; i<cardBaseV; ++i) {
            res += (U(2*i) - V(2*i))*(U(2*i) - V(2*i)) + (U(2*i+1) - V(2*i+1))*(U(2*i+1) - V(2*i+1));
        }
        return res;
    }

} // namespace procedure

#endif