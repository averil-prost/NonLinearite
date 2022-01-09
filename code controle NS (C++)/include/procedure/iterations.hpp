#ifndef __ITERATIONS__
#define __ITERATIONS__

#include "../inclusions.hpp"
#include "../datastructure.hpp"
#include "../builder.hpp"
#include "../interface.hpp"
#include "../visu.hpp"

namespace procedure {

    // ajoute le vecteur ctilde (XX,YY,.) à vect
    void addbicontract (Vecteur& vect, const datastructure::Matrice3D& ctilde, const Vecteur& XX, const Vecteur& YY, const datastructure::bordermgr& bman, bool verbose) {
        if (verbose) std::cout << "Entrée dans add bicontract" << std::endl;
        
        //   La coordonnée l de la partie m (1 ou 2) du vecteur res est donnée par
        //   
        //       somme (i=1,2) somme (j=1,n) somme (k=1,n) <partie i de u>_j * <partie m de v>_k * c_i (j,k,l)
        //   
        //   où  - partie alpha d'un vecteur : alpha=1, les x, alpha=2, les y, d'abord x puis y 
        //       - u_j : je coordonnée de u dans la base associée
        //       - c_i : ctilde.first si i=1 (dérivées en x), ctilde.second si i=2 (dérivées en y)
        
        int indk=0;
        int ielem1=0, nelem1=0, col1=0;
        int ielem2=0, nelem2=0, col2=0;
        double linevalue=0.0, cdot_YY1=0.0, cdot_YY2=0.0;
        for (auto kvalue : ctilde.data) { // pour chaque coordonnée k
            if (bman.isInside.at(indk)) { // la fonction test est dans H^1_0
                cdot_YY1=0.0; 
                cdot_YY2=0.0;
                for (int line : ctilde.nonzeros.at(indk)) { // faire le produit kvalue.first * <partie 1 de XX> + kvalue.second * <partie 2 de XX> 
                    nelem1 = kvalue.first.getNbTermL (line);
                    nelem2 = kvalue.second.getNbTermL (line);
                    linevalue = 0.0; ielem1=0; ielem2=0;
                    while ((ielem1 < nelem1) || (ielem2 < nelem2)) {
                        col1 = kvalue.first.getIndLC (line, ielem1);
                        col2 = kvalue.second.getIndLC (line, ielem2);
                        if (col1==col2) {
                            linevalue += kvalue.first (line,col1) * XX(2*col1) + kvalue.second(line,col2) * XX(2*col2+1);
                            ++ielem1; ++ielem2;
                        } else if (col1 < col2) {
                            linevalue += kvalue.first (line,col1) * XX(2*col1);
                            ++ielem1;
                        } else {
                            linevalue += kvalue.second(line,col2) * XX(2*col2+1);
                            ++ielem2;
                        }
                    }
                    cdot_YY1 += linevalue * YY(2*line);        // faire le produit scalaire <produit> * <partie 1 de YY>
                    cdot_YY2 += linevalue * YY(2*line+1);      // faire le produit scalaire <produit> * <partie 2 de YY>
                }
                vect.set (2*indk,   vect(2*indk  ) + cdot_YY1);
                vect.set (2*indk+1, vect(2*indk+1) + cdot_YY2);
            }
            ++indk;
        }
        if (verbose) std::cout << "Sortie de add bicontract" << std::endl;
    }

    // DOCUMENTATION TIME les différentes contractions
    // Notons w ce qu'on veut calculer, u ce qu'on a, v la fonction test.
    // Dans les schémas apparaissent ctilde (w, u, v), ctilde (v, u, w) et ctilde (u, v, w).
    // La position de v donne l'indice des lignes, celle de w l'indice des colonnes.
    // On doit donc ressortir respectivement C_{k,i}, C_{i,k} et C_{j,k}, où l'indice non présent
    // est l'indice de somme (exemple : C_{k,j} := (sum_i=1^n ctilde (i,j,k))_{k,j}).

    // convention de nomenclature : les noms i,j,k font référence aux indices d'accès à ctilde (i,j,k)
    // le premier indice donne les lignes de la matrice contractée, le deuxième les colonnes.

    // dernière remarque : dans l'implémentation, on a toujours line=i, col=j, indk=k.
    // ah, et la pseudo-élimination concerne les *lignes* de la sortie, donc ça change à chaque fois.

    // ajoute la matrice C_{k,i} := (sum_j=1^n u_j ctilde (i,j,k))_{k,i} à mat
    void addcontract_ki (MatriceL& mat, const datastructure::Matrice3D& ctilde, const Vecteur& XX, const datastructure::bordermgr& bman, bool verbose) {
        if (verbose) std::cout << "Entrée dans addcontract_ki" << std::endl;
        int indk=0; // parcourt les lignes
        int ielem1=0, nelem1=0, col1=0;
        int ielem2=0, nelem2=0, col2=0;
        double value = 0.0;
        for (auto kvalue : ctilde.data) { // paire des matrices ctilde_dx et ctilde_dy
            if (bman.isInside.at(indk)) { // ne pas toucher à la pseudo-élimination
                for (int line : ctilde.nonzeros.at(indk)) { // parcourt les i
                    nelem1 = kvalue.first.getNbTermL (line);
                    nelem2 = kvalue.second.getNbTermL (line);
                    value = 0.0; ielem1=0; ielem2=0;
                    while ((ielem1 < nelem1) || (ielem2 < nelem2)) {
                        col1 = kvalue.first.getIndLC (line, ielem1);
                        col2 = kvalue.second.getIndLC (line, ielem2);
                        if (col1==col2) {
                            value += kvalue.first (line,col1) * XX(2*col1) + kvalue.second(line,col2) * XX(2*col2+1);
                            ++ielem1; ++ielem2;
                        } else if (col1 < col2) {
                            value += kvalue.first (line,col1) * XX(2*col1);
                            ++ielem1;
                        } else {
                            value += kvalue.second(line,col2) * XX(2*col2+1);
                            ++ielem2;
                        }
                    }
                    mat.set (2*indk,   2*line,   mat(2*indk,   2*line  ) + value);
                    mat.set (2*indk+1, 2*line+1, mat(2*indk+1, 2*line+1) + value);
                } // for line
            } // if k isInside
            ++indk;
        } // for kvalue
        if (verbose) std::cout << "Sortie de addcontract_ki" << std::endl;
    }

    // ajoute la matrice C_{i,k} := (sum_j=1^n u_j ctilde (i,j,k))_{i,k} à mat
    void addcontract_ik (MatriceL& mat, const datastructure::Matrice3D& ctilde, const Vecteur& XX, const datastructure::bordermgr& bman, bool verbose) {
        if (verbose) std::cout << "Entrée dans addcontract_ik" << std::endl;
        int indk=0; // parcourt les lignes
        int ielem1=0, nelem1=0, col1=0;
        int ielem2=0, nelem2=0, col2=0;
        double value = 0.0;
        for (auto kvalue : ctilde.data) { // paire des matrices ctilde_dx et ctilde_dy
            for (int line : ctilde.nonzeros.at(indk)) { // parcourt les colonnes
                if (bman.isInside.at(line)) { // ne pas toucher à la pseudo-élimination
                    nelem1 = kvalue.first.getNbTermL (line);
                    nelem2 = kvalue.second.getNbTermL (line);
                    value = 0.0; ielem1=0; ielem2=0;
                    while ((ielem1 < nelem1) || (ielem2 < nelem2)) {
                        col1 = kvalue.first.getIndLC (line, ielem1);
                        col2 = kvalue.second.getIndLC (line, ielem2);
                        if (col1==col2) {
                            value += kvalue.first (line,col1) * XX(2*col1) + kvalue.second(line,col2) * XX(2*col2+1);
                            ++ielem1; ++ielem2;
                        } else if (col1 < col2) {
                            value += kvalue.first (line,col1) * XX(2*col1);
                            ++ielem1;
                        } else {
                            value += kvalue.second(line,col2) * XX(2*col2+1);
                            ++ielem2;
                        }
                    }
                    mat.set (2*line,   2*indk,   mat(2*line,   2*indk  ) + value);
                    mat.set (2*line+1, 2*indk+1, mat(2*line+1, 2*indk+1) + value);
                } // if i isInside
            } // for line
            ++indk;
        } // for kvalue
        if (verbose) std::cout << "Sortie de addcontract_ik" << std::endl;
    }

    // ajoute la matrice C_{j,k} := (sum_i=1^n u_i ctilde (i,j,k))_{j,k} à mat
    void addcontract_jk (MatriceL& mat, const datastructure::Matrice3D& ctilde, const Vecteur& XX, const datastructure::bordermgr& bman, bool verbose) {
        if (verbose) std::cout << "Entrée dans addcontract_jk" << std::endl;
        int indk=0;
        int ielem=0, nelem=0, col=0;
        double XXvalue = 0.0;
        Vecteur value (ctilde.cubeSize);
        // le fun de cette permutation, c'est que le stockage sparse est dans le mauvais sens
        for (auto kvalue : ctilde.data) { // paire des matrices ctilde_dx et ctilde_dy
            for (int line : ctilde.nonzeros.at(indk)) { // parcourt les i
                // somme sur la première dimension => dérivée par rapport à x
                XXvalue = XX (2*line);
                nelem = kvalue.first.getNbTermL (line);
                for (ielem=0; ielem<nelem; ++ielem) {
                    col = kvalue.first.getIndLC (line, ielem);
                    value.set (col, value(col) + XXvalue * kvalue.first (line, col)); 
                }
                // somme sur la deuxième dimension => dérivée par rapport à y
                XXvalue = XX (2*line+1);
                nelem = kvalue.second.getNbTermL (line);
                for (ielem=0; ielem<nelem; ++ielem) {
                    col = kvalue.second.getIndLC (line, ielem);
                    value.set (col, value(col) + XXvalue * kvalue.second (line, col)); 
                }
            } // for line
            // recopie des valeurs pour chaque j sauf si pseudo-élimination
            for (col=0; col<ctilde.cubeSize; ++col) { 
                if (bman.isInside.at(col)) {  
                    mat.set (2*col,   2*indk,   mat(2*col,   2*indk  ) + value(col));
                    mat.set (2*col+1, 2*indk+1, mat(2*col+1, 2*indk+1) + value(col));
                }
                // DANGEREUX mais opti : reset de value en même temps
                value.set(col,0.0);
            } // for j if j isInside
            ++indk;
        } // for kvalue
        if (verbose) std::cout << "Sortie de addcontract_jk" << std::endl;
    }

    // itérations du problème direct
    Vecteur iterate_forward (datastructure::insta& insta, Vecteur& Un, Vecteur& fn, Vecteur& dirichnp1, Vecteur& dirichn,
        const datastructure::bordermgr& bman, double deltaT, double nu, bool checkStationnary, bool verbose) {
        if (verbose) std::cout << "Entrée dans iterate forward" << std::endl;

        MatriceL *tobeinverted; 
        MatriceL current_tbi (insta.tobeinverted.getDimX(), insta.tobeinverted.getDimY());
        Vecteur scdmembre (2*insta.cardBaseV+insta.cardBaseP+1);

        /////////////////////////////////////////////////////////
        // calcul du second membre
        /////////////////////////////////////////////////////////
        int col=0,ielem=0,nelem=0;
        double linevalue1=0.0, linevalue2=0.0, matvalue=0.0;
        for (int line : bman.memory_to_mesh) { // ne mettre à jour que les noeuds intérieurs
            linevalue1 = 0.0; linevalue2=0.0;
            // contribution des termes multipliés par la masse
            nelem = insta.uv.getNbTermL (line);
            for (ielem=0; ielem<nelem; ++ielem) {
                col = insta.uv.getIndLC (line, ielem);
                matvalue = insta.uv (line, col);
                linevalue1 += matvalue*((dirichnp1(2*col  )+Un(2*col  )-dirichn(2*col  )) / deltaT + fn (2*col));
                linevalue2 += matvalue*((dirichnp1(2*col+1)+Un(2*col+1)-dirichn(2*col+1)) / deltaT + fn (2*col+1));
            }
            // contribution des termes multipliés par la rigidité
            nelem = insta.gradugradv.getNbTermL (line);
            for (ielem=0; ielem<nelem; ++ielem) {
                col = insta.gradugradv.getIndLC (line, ielem);
                matvalue = nu * insta.gradugradv(line, col);
                linevalue1 -= matvalue*dirichnp1(2*col);
                linevalue2 -= matvalue*dirichnp1(2*col+1);
            }
            scdmembre.set (2*line  , linevalue1);
            scdmembre.set (2*line+1, linevalue2);
        }
        // contribution du relèvement à la pression
        // !!!!!! on dispose de mdivu1u2, on veut sa transposée => sparse dans le sens des colonnes
        // line parcourt les lignes de dirich
        for (int line=0; line<2*insta.cardBaseV; ++line) {
            nelem = insta.mdivu1v2.getNbTermL (line);
            linevalue1 = dirichnp1 (line);
            for (ielem=0; ielem<nelem; ++ielem) {
                col = insta.mdivu1v2.getIndLC (line,ielem);
                scdmembre.set(2*insta.cardBaseV+col, scdmembre(2*insta.cardBaseV+col) - linevalue1 * insta.mdivu1v2 (line, col));
            }
        }

        /////////////////////////////////////////////////////////
        // calcul de la matrice à résoudreLeSystèmeLinéaire-er
        /////////////////////////////////////////////////////////
        if (insta.name == "navierim") {
            // current_tbi = insta.tobeinverted;
            builder::copyMatriceL (insta.tobeinverted, current_tbi);
            procedure::addcontract_ki (current_tbi, insta.ctilde, Un+dirichn, bman, verbose);
            procedure::addbicontract (scdmembre, insta.ctilde, (Un+dirichn)*(-1.0), dirichnp1, bman, verbose);
            tobeinverted = &current_tbi;
            // tobeinverted = &insta.tobeinverted;
        }
        else if (insta.name == "navierex") {
            procedure::addbicontract (scdmembre, insta.ctilde, (Un+dirichn)*(-1.0), Un+dirichn, bman, verbose);
            tobeinverted = &insta.tobeinverted;
        }
        else { // DEFAULT : STOKES
            tobeinverted = &insta.tobeinverted;
        }

        // TRÈS UTILE si la solution se stabilise : sinon, NaN direct
        double normeInit=0.0; if (checkStationnary) normeInit = ((*tobeinverted) * Un - scdmembre).norme(); 

        // résolution du système linéaire
        if (!checkStationnary || normeInit > 1e-8) {
            if (insta.name != "stokes") {
                return (*tobeinverted).solveIte (scdmembre, Un, "GMRES", "Id", 1.47, 8000, 500, 1e-8, 0);
            } else {
                return (*tobeinverted).solveIte (scdmembre, Un, "KrylovSym", "Id", 1.47, 8000, 8000, 1e-8, 0);
            }
        }
        // if stationary solution reached
        return Un;
    }

    // itérations du problème adjoint
    // un seul schéma, pas de relèvement (toujours nul au bord, quand même une pseudo-élimination)
    Vecteur iterate_backward (datastructure::insta& insta, Vecteur& Wnp1, Vecteur& Unp1, Vecteur& U_Qnp1,
        const datastructure::bordermgr& bman, double deltaT, bool verbose) {
        if (verbose) std::cout << "Entrée dans iterate backward" << std::endl;

        Vecteur scdmembre (2*insta.cardBaseV+insta.cardBaseP+1);

        /////////////////////////////////////////////////////////
        // calcul du second membre
        /////////////////////////////////////////////////////////
        int col=0,ielem=0,nelem=0;
        double linevalue1=0.0, linevalue2=0.0, matvalue=0.0;
        for (int line : bman.memory_to_mesh) { // ne mettre à jour que les noeuds intérieurs
            linevalue1=0.0; linevalue2=0.0;
            // contribution des termes multipliés par la masse
            nelem = insta.uv.getNbTermL (line);
            for (ielem=0; ielem<nelem; ++ielem) {
                col = insta.uv.getIndLC (line, ielem);
                matvalue = insta.uv (line, col);
                linevalue1 += matvalue*(Wnp1(2*col  ) / deltaT + Unp1 (2*col  ) - U_Qnp1(2*col  ));
                linevalue2 += matvalue*(Wnp1(2*col+1) / deltaT + Unp1 (2*col+1) - U_Qnp1(2*col+1));
            }
            scdmembre.set (2*line  , linevalue1);
            scdmembre.set (2*line+1, linevalue2);
        }

        /////////////////////////////////////////////////////////
        // calcul de la matrice à résoudreLeSystèmeLinéaire-er
        /////////////////////////////////////////////////////////
        if (insta.name != "stokes") {
            MatriceL current_tbi (insta.tobeinverted.getDimX(), insta.tobeinverted.getDimY());
            builder::copyMatriceL (insta.tobeinverted, current_tbi);
            // current_tbi = insta.tobeinverted;
            // pour le futur : faire une horreur qui combine les deux ?
            procedure::addcontract_ik (current_tbi, insta.ctilde, Unp1, bman, verbose);
            procedure::addcontract_jk (current_tbi, insta.ctilde, Unp1, bman, verbose);

            // std::cout << "\t\t\tnorme courante : " << (current_tbi * Wnp1 - scdmembre).norme () << std::endl;

            return current_tbi.solveIte (scdmembre, Wnp1, "GMRES", "Id", 1.47, 8000, 500, 1e-8, 0);
        }

        // IF STOKES (moche, mais évite une matrice temporaire)
        // return insta.tobeinverted.solveIte (scdmembre, Wnp1, "GMRES", "Id", 1.47, 8000, 500, 1e-8, 0);
        return insta.tobeinverted.solveIte (scdmembre, Wnp1, "KrylovSym", "Id", 1.47, 8000, 8000, 1e-7, 0);
    }

} // namespace procedure

#endif