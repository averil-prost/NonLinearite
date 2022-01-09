#ifndef __DATASTRUCTURE__
#define __DATASTRUCTURE__

#include "inclusions.hpp"

namespace datastructure {
    
    struct Matrice3D {
        // données
        int cubeSize;
        // toutes vont être non vides. 
        // Première contient les dérivées par rapport à x, 2e par rapport à y
        std::vector<std::pair<MatriceL,MatriceL>> data; // indice : k, le 3eme
		std::vector<std::vector<int>> nonzeros; // pour chaque k, liste des indices des j non vides
    };

    struct extension_ef_d2 {
        element_fini_D2 *ef_v, *ef_p;
        int cardBaseVref, cardBasePref;
        // deux composantes de c (qui sera évaluée en ctilde)
        double *** cDx;
        double *** cDy;
        // deux composantes des matrices de croisement entre vitesse et pression 
        double ** bDx; // base vitesse en ligne, base pression en colonne
        double ** bDy;
        // intégrales des fonctions de base de la pression
        double * intp;
    };

    // pour les schémas instationnaires
    struct insta {
        std::string name; // stokes, navierex, navierim 
        MatriceL tobeinverted; // partie stable
        MatriceL uv;
        MatriceL gradugradv;
        Matrice3D ctilde;
        MatriceL mdivu1v2;

        int cardBaseV, cardBaseP;
    };

    // schéma stokes stationnaire : pour condition initiale
    struct sta {
        MatriceL tobeinverted;
    };

    struct schemas {
        datastructure::insta insta;
        datastructure::sta sta;
    };

    struct bordermgr {
        std::vector<bool> isInside; // vecteur de booléens : true si à l'intérieur, false si au bord
        int nBorderNodes;
        std::vector<int> mesh_to_memory; // numérotation du maillage à l'ordre en mémoire (sans les bords)
        std::vector<int> memory_to_mesh; // de l'ordre (contigue) en mémoire à la numérotation du maillage
        // exemple : si le maillage a pour noeuds {1, 2, 3, 4}, et que 3 est au bord
        // mesh_to_memory sera {1, 2, -1, 3}, memory_to_mesh {1, 2, 4}
        // d'où le parcours de memory_to_mesh ira directement itérer sur les noeuds de l'intérieur.

        bordermgr (int nnodes) 
            : isInside (nnodes, true)
            {}
    };

} // namespace datastructure

#endif