#ifndef __BUILDER__
#define __BUILDER__

#include "inclusions.hpp"
#include "datastructure.hpp"

namespace builder {

    // encapsule (aaaa) toutes les infos nécessaires à la gestion des bords
    datastructure::bordermgr bordermgr (Meshes& mesh, int cardBaseVRef) {
        datastructure::bordermgr bman (mesh.getNNodes()); // initialisation de isInside à true
        std::vector<Vecteur*> nodes = mesh.getNodes();
        bman.nBorderNodes=0;
        // détection des bords dans le maillage grossier (juste triangles)
        for (auto arc : mesh.getMesh1D ()) { // si c'est là-dedans, c'est à la frontière
                if (bman.isInside.at(arc[1])) bman.nBorderNodes++; // si on considérait encore qu'il était true
                if (bman.isInside.at(arc[2])) bman.nBorderNodes++;
                bman.isInside.at(arc[1]) = false; // frontière pas forcément orientée... 
                bman.isInside.at(arc[2]) = false; // obligée de faire 2 tests !
        }
        // ajout des points dus au maillage d'ordre > 1 
        int inode=0, ibegin=0, iend=0;
        double deltax=0, deltay=0, xdep=0, ydep=0;
        for (auto elem : mesh.getMesh2D()) { // pour chaque triangle (pas fou, mais quoi d'autre ? ...)
            for (ibegin=1, iend=2; ibegin<=3; ++ibegin, iend = iend%3 + 1) { // permutation [1,2], [2,3], [3,1] (hihi)
                if ((!bman.isInside.at(elem[ibegin])) && (!bman.isInside.at(elem[iend]))) { // arête [ibegin, iend]
                    xdep   = (*nodes.at(elem[ibegin]))(0);          ydep   = (*nodes.at(elem[ibegin]))(1);
                    deltax = (*nodes.at(elem[iend]))(0) - xdep;     deltay = (*nodes.at(elem[iend]))(1) - ydep;
                    for (inode=4; inode<cardBaseVRef+1; ++inode) { // pour tout autre point de l'élément que les coins
                        // test de colinéarité : det de la matrice [[x - xdep, deltax], [y - ydep, deltay]] ~ 0
                        if (std::fabs(((*nodes.at(elem[inode]))(0) - xdep)*deltay - ((*nodes.at(elem[inode]))(1) - ydep)*deltax) < 1e-10) {
                            // enregistrement du noeud en question comme étant sur la frontière
                            if (bman.isInside.at(elem[inode])) bman.nBorderNodes++;
                            bman.isInside.at(elem[inode]) = false; 
                        } 
                    } // for potentiels noeuds sur l'arête
                } // if les deux extrémités sont sur la frontière
            } // for chaque arête du triangle
        } // for chaque triangle
        // écriture des conversions
        int innernode=0, node=0;
        for (bool is : bman.isInside) {
            if (is) { // noeud à l'intérieur
                bman.memory_to_mesh.push_back (node);       // indice du noeud selon gmsh
                bman.mesh_to_memory.push_back (innernode);  // dernier indice (mémoire) libre
                innernode++;
            } else { // noeud sur le bord : pas d'image par mesh_to_memory, mais faut mettre qqch (redondant avec isInside, je sais)
                bman.mesh_to_memory.push_back (-1);
            }
            node++;
        }
        return bman;
    }
    
    // renvoie le triangle de Pascal jusqu'à dim
    // k parmi n (si n<=dim) : res.at(n).at(k)
    std::vector<std::vector<int>> triPascal (const int dim) {
        std::vector<std::vector<int>> res;
        int k=0, n=0;
        for (n=0; n<=dim; ++n) {
            std::vector<int> line;
            line.push_back(1);
            for (k=1; k<n; ++k) {
                line.push_back ((*(res.end()-1)).at(k-1) + (*(res.end()-1)).at(k));
            }
            if (n>0) line.push_back(1);
            res.push_back(line);
        }
        return res;
    }

    // renvoie les intégrales sur le triangle de ref des x^a*y^b
    // ordre : for a=0 to degMax for b=0 to degMax-a
    std::vector<double> getIntegralesMonomes (int degMax) {
        std::vector<double> res;
        int nbPoints = (degMax+1)*(degMax+2)/2, a=0, b=0, i=0;
        double value=0.0;
        // compute pascal triangle up to max degree + 1
        std::vector<std::vector<int>> tripasc = builder::triPascal (degMax+1);
        for (a=0; a<=degMax; ++a) { 
            for (b=0; b<=degMax - a; ++b) { // fonction x^a * y^b
                value=0;
                for (i=0; i<=b+1; ++i) {
                    value += tripasc.at(b+1).at(i) *(i%2==0?1.0:-1.0) / (i + 1 + a);
                }
                value *= 1.0/(b+1);
                res.push_back (value);
            }
        }
        return res;
    }

    // ajoute les matrices élémentaires propres à ce code 
    datastructure::extension_ef_d2 extension_ef_d2 (element_fini_D2* ef_v, element_fini_D2* ef_p) {
        datastructure::extension_ef_d2 extef;
        extef.ef_v = ef_v;
        extef.ef_p = ef_p;
        extef.cardBaseVref = extef.ef_v->getTaille(); 
        extef.cardBasePref = extef.ef_p->getTaille(); 
        int degMaxAIntegrer=3*extef.ef_v->getOrdre(); // HYPOTHÈSE le cardinal de la base pression est inférieur à celui de la base vitesse 
        int nbCoeffV=(extef.ef_v->getOrdre()+1)*(extef.ef_v->getOrdre()+2)/2;
        int nbCoeffP=(extef.ef_p->getOrdre()+1)*(extef.ef_p->getOrdre()+2)/2;
        int i=0,j=0,k=0,coeffi=0,coeffj=0,coeffk=0,degTotalX=0,degTotalY=0;
        double produitCoeffs=0.0;
        poly_2D *phi_i, *phi_j, *phi_k;

        // récupération des intégrales des x^a*y^b pour a et b <= degMaxAIntegrer
        std::vector<double> integralesMonomes = builder::getIntegralesMonomes (degMaxAIntegrer);        
        std::vector<std::vector<int>> corresp; // corresp.at(a).at(b) donne l'indice de intégrale(x^a*y^b)
        k=0;
        for (i=0; i<=degMaxAIntegrer; ++i) {
            std::vector<int> line;
            for (j=0; j<=degMaxAIntegrer-i; ++j) {
                line.push_back(k);
                ++k;
            }
            corresp.push_back(line);
        }

        // évaluation des matrices de référence pour ctilde
        extef.cDx = new double**[extef.cardBaseVref];
        extef.cDy = new double**[extef.cardBaseVref];
        for (i=0; i<extef.cardBaseVref; ++i) {
            phi_i = extef.ef_v->getFonctLagrange(i);
            extef.cDx[i] = new double*[extef.cardBaseVref];
            extef.cDy[i] = new double*[extef.cardBaseVref];
            for (j=0; j<extef.cardBaseVref; ++j) {
                phi_j = extef.ef_v->getFonctLagrange(j);
                extef.cDx[i][j] = new double[extef.cardBaseVref-i];
                extef.cDy[i][j] = new double[extef.cardBaseVref-i];
                for (k=i; k<extef.cardBaseVref; ++k) {
                    phi_k = extef.ef_v->getFonctLagrange (k);
                    extef.cDx[i][j][k-i]=0.0;
                    extef.cDy[i][j][k-i]=0.0;

                    for (coeffi=0; coeffi<nbCoeffV; ++coeffi) {
                        for (coeffj=0; coeffj<nbCoeffV; ++coeffj) {
                            for (coeffk=0; coeffk<nbCoeffV; ++coeffk) {
                                produitCoeffs = phi_i->getMonome2D(coeffi)->getCoeff() * phi_j->getMonome2D(coeffj)->getCoeff() * phi_k->getMonome2D(coeffk)->getCoeff();
                                degTotalX = phi_i->getMonome2D(coeffi)->getPowX() + phi_j->getMonome2D(coeffj)->getPowX() + phi_k->getMonome2D(coeffk)->getPowX();
                                degTotalY = phi_i->getMonome2D(coeffi)->getPowY() + phi_j->getMonome2D(coeffj)->getPowY() + phi_k->getMonome2D(coeffk)->getPowY();
                                
                                if (degTotalX > 0)
                                    extef.cDx[i][j][k-i] += produitCoeffs * phi_j->getMonome2D(coeffj)->getPowX() * integralesMonomes.at (corresp.at (degTotalX - 1).at (degTotalY));
                                if (degTotalY > 0)
                                    extef.cDy[i][j][k-i] += produitCoeffs * phi_j->getMonome2D(coeffj)->getPowY() * integralesMonomes.at (corresp.at (degTotalX).at (degTotalY - 1));
                            }
                        }
                    }
                }
            }
        }

        // ORDRE DU TREILLIS

        // évaluation des matrices de référence pour b
        extef.bDx = new double*[extef.cardBaseVref];
        extef.bDy = new double*[extef.cardBaseVref];
        for (i=0; i<extef.cardBaseVref; ++i) {
            phi_i = extef.ef_v->getFonctLagrange(i);
            extef.bDx[i] = new double[extef.cardBasePref];
            extef.bDy[i] = new double[extef.cardBasePref];
            for (j=0; j<extef.cardBasePref; ++j) {
                phi_j = extef.ef_p->getFonctLagrange(j);
                extef.bDx[i][j]=0.0;
                extef.bDy[i][j]=0.0;

                for (coeffi=0; coeffi<nbCoeffV; ++coeffi) {
                    for (coeffj=0; coeffj<nbCoeffP; ++coeffj) {
                        produitCoeffs = phi_i->getMonome2D(coeffi)->getCoeff() * phi_j->getMonome2D(coeffj)->getCoeff();
                        degTotalX = phi_i->getMonome2D(coeffi)->getPowX() + phi_j->getMonome2D(coeffj)->getPowX();
                        degTotalY = phi_i->getMonome2D(coeffi)->getPowY() + phi_j->getMonome2D(coeffj)->getPowY();
                        
                        if (degTotalX > 0)
                            extef.bDx[i][j] -= produitCoeffs * phi_i->getMonome2D(coeffi)->getPowX() * integralesMonomes.at (corresp.at (degTotalX - 1).at (degTotalY));
                        if (degTotalY > 0)
                            extef.bDy[i][j] -= produitCoeffs * phi_i->getMonome2D(coeffi)->getPowY() * integralesMonomes.at (corresp.at (degTotalX).at (degTotalY - 1));
                    }
                }
            }
        }

        // calcul des intégrales des fonctions de base de la pression
        extef.intp = new double [extef.cardBasePref];
        for (j=0; j<extef.cardBasePref; ++j) {
            phi_j = extef.ef_p->getFonctLagrange(j);
            extef.intp [j] = 0.0;
            for (coeffj=0; coeffj<nbCoeffP; ++coeffj) {
                extef.intp [j] += phi_j->getMonome2D(coeffj)->getCoeff() * \
                    integralesMonomes.at(corresp.at(phi_j->getMonome2D(coeffj)->getPowX()).at(phi_j->getMonome2D(coeffj)->getPowY()));
            }
        }

        return extef;
    } 

    void interpolate (Vecteur& interp, Vecteur source (Vecteur&, double), Meshes& mesh_v, double t) {
        std::cout << "Entrée dans interpolate" << std::endl;
        int inode=0;
        Vecteur value (3);
        std::vector<Vecteur*> nodes = mesh_v.getNodes();
        for (auto node : nodes) { // tous les noeuds du maillage vitesse (ordre >= 1)
            value = source (*node, t);
            interp.set(2*inode,   value(0));
            interp.set(2*inode+1, value(1));
            ++inode;
        }
        std::cout << "Sortie de interpolate" << std::endl;
    }

    // loc = 0 si n est pair, 1 si n impair
    void relever (std::vector<Vecteur>& dirich, Vecteur dirichlet(Vecteur&, double), Meshes& mesh_v, 
            int loc, double t, int cardBaseV, const datastructure::bordermgr& bman) {
        Vecteur value (3);
        for (int k=0; k<cardBaseV; ++k) {
            if (!bman.isInside.at(k)) {
                value = dirichlet (*mesh_v.getNodes().at(k), t);
                dirich.at(loc).set(2*k, value(0)); 
                dirich.at(loc).set(2*k+1, value(1)); 
            }
        }
    }

    // l'ordre des fonctions de base dans les données du maillage n'est pas forcément évident.
    // on a donc besoin de correspondance : donne la permutation de {1,...,n} telle que 
    // le numéro de la fonction i soit correspondance[i], où i parcourt l'ordre donné par le treillis

    // correspondance : ordre_treillis -> ordre_maillage (de gmsh)
    std::vector<int> get_corresp (element_fini_D2& ef, Meshes& mesh) {
        std::vector<int> corresp; 
        std::vector<double> distances; std::vector<double>::iterator mindist;

        double ** treillis = ef.getTreilli(); // treillis des noeuds des fonctions de base
        int itaille=0, taille=ef.getTaille(); // nombre de fonctions de base sur un triangle de référence
        std::vector<Vecteur*> noeudsPremierTriangle;
        // récupérer la transformation du premier triangle
        Vecteur imageDeRef (3), p1 (3), p2mp1 (3), p3mp1 (3);
        p1 = *mesh.getNodes()[mesh.getMesh2D()[0][1]];
        p2mp1 = *mesh.getNodes()[mesh.getMesh2D()[0][2]] - p1; 
        p3mp1 = *mesh.getNodes()[mesh.getMesh2D()[0][3]] - p1; 
   
        for (itaille=0; itaille<taille; ++itaille) { // récupérer les noeuds des fonctions de base associées au premier élément
            noeudsPremierTriangle.push_back(mesh.getNodes()[mesh.getMesh2D()[0][itaille+1]]);
        } 

        for (itaille=0; itaille<taille; ++itaille) { // pour chaque fonction de base de référence dans l'ordre du treillis
            imageDeRef = p1 + p2mp1 * treillis[itaille][0] + p3mp1 * treillis[itaille][1]; // récupérer l'image du noeud de ref
            for (Vecteur* pgmshfunc : noeudsPremierTriangle) { // calculer les distances aux points de gmsh
                distances.push_back ((imageDeRef - *pgmshfunc).norme());
            }
            mindist = std::min_element (distances.begin(), distances.end()); // récupérer la plus petite distance
            if (*mindist > 1e-8) {std::cout << "ATTENTION pas de correspondance pour le point " << itaille << ", mindist=" << *mindist << std::endl;}
            else {
                corresp.push_back (std::distance(distances.begin(), mindist)); // indice du minimum
            }
            distances.clear();
        }
        return corresp;
    }

    // renvoie les coefficients de la projection de f(t) sur l'espace associé au maillage
    Vecteur project_f (double f (Vecteur&, double), double t, Meshes& mesh, element_fini_D2& ef, const std::vector<int>& correspondance) {
        // masse des éléments finis considérés : mesh est déjà converti
        Bilinear myBilinear;
        MatriceL masse (mesh.getNNodes());
        myBilinear.uv2D (mesh.getNodes(), mesh.getMesh2D(), ef, &masse);
        // projection
        int nnodes = mesh.getNNodes(), nquad=ef.getNbPointQuad(), taille=ef.getTaille(), iquad=0,ifunc=0,itaille=0;
        double** quad = ef.getPointQuad(); // récupérer la super quadrature
        double factor=0.0, jacobien=0.0; // poids * f(point de quadrature dans le triangle d'origine) * |jacobien| de F : K^ -> K
        std::vector<Vecteur*> nodes = mesh.getNodes();
        Vecteur dprods (nnodes), currentpoint(3);
        for (int* elem : mesh.getMesh2D()) { // pour chaque triangle du maillage
            jacobien = std::fabs(((*nodes[elem[2]])(0)-(*nodes[elem[1]])(0))*((*nodes[elem[3]])(1)-(*nodes[elem[1]])(1)) \
                               - ((*nodes[elem[2]])(1)-(*nodes[elem[1]])(1))*((*nodes[elem[3]])(0)-(*nodes[elem[1]])(0)));
            for (iquad=0; iquad<nquad; ++iquad) { // pour chaque noeud de la formule de quadrature
                // calcul du point courant
                currentpoint = *nodes[elem[1]] + (*nodes[elem[2]] - *nodes[elem[1]])*quad[iquad][0] \
                                               + (*nodes[elem[3]] - *nodes[elem[1]])*quad[iquad][1];
                // point de quadrature * f(x) * jacobien
                factor = quad[iquad][2] * f(currentpoint,t) * jacobien;
                for (itaille=0; itaille<taille; ++itaille) { // pour toutes les fonctions chapeau non nulles sur le triangle
                    // getValLagrange prend le numéro de la fonction dans {1,taille}
                    dprods.set(elem[correspondance[itaille]+1], dprods(elem[correspondance[itaille]+1]) + factor * ef.getValLagrange(itaille,iquad));
                }
            }
        }
        // std::cout << "Produits scalaires : " << visu::tostring (dprods);
        return masse.solveIte (dprods, dprods*0.0, (char*)"CG", (char*)"ILU", 1.47, 3000, 3000, 1e-10, 0);
    }

    // version projection
    Vecteur getSol_proj (double solv1 (Vecteur&, double), double solv2 (Vecteur&, double), double solp (Vecteur&, double), double t, 
        Meshes& mesh_v, Meshes& mesh_p, element_fini_D2& ef_v,  element_fini_D2& ef_p, int cardBaseV, int cardBaseP) {
        // plus de calcul mais moins de dépendances cachées dans le code...
        std::vector<int> correspondance_v = builder::get_corresp (ef_v, mesh_v); 
        std::vector<int> correspondance_p = builder::get_corresp (ef_p, mesh_p);
        Vecteur sol (2*cardBaseV + cardBaseP + 1), projv (cardBaseV), projp (cardBaseP);
        int i=0;
        // coordonnée x
        projv = builder::project_f (solv1, t, mesh_v, ef_v, correspondance_v);
        for (i=0; i<cardBaseV; ++i) {
            sol.set (2*i, projv (i));
        }
        // coordonnée y
        projv = builder::project_f (solv2, t, mesh_v, ef_v, correspondance_v);
        for (i=0; i<cardBaseV; ++i) {
            sol.set (2*i+1, projv (i));
        }
        // pression
        projp = builder::project_f (solp, t, mesh_p, ef_p, correspondance_p);
        for (i=0; i<cardBaseP; ++i) {
            sol.set (2*cardBaseV+i, projp (i));
        }
        return sol;
    }

    // version interpolation
    Vecteur getSol_interp (double solv1 (Vecteur&, double), double solv2 (Vecteur&, double), double solp (Vecteur&, double), double t, 
        Meshes& mesh_v, Meshes& mesh_p, int cardBaseV, int cardBaseP) {

        Vecteur sol (2*cardBaseV + cardBaseP + 1);
        int inode=0; 
        for (Vecteur* node : mesh_v.getNodes()) { // pour chaque noeud du maillage
            sol.set(2*inode  , solv1(*node, t));
            sol.set(2*inode+1, solv2(*node, t));
            ++inode;
        }

        inode=0;
        for (Vecteur* node : mesh_p.getNodes()) { // pour chaque noeud du maillage
            sol.set(2*cardBaseV+inode, solp(*node, t));
            ++inode;
        }
        return sol;
    }

    void free (datastructure::extension_ef_d2& extef) {
        int i=0, j=0;
        for (i=0; i<extef.cardBaseVref; ++i) {
            for (j=0; j<extef.cardBaseVref; ++j) {
                std::free (extef.cDx[i][j]);
                std::free (extef.cDy[i][j]);
            }
            std::free (extef.cDx[i]);
            std::free (extef.cDy[i]);

            std::free (extef.bDx[i]);
            std::free (extef.bDy[i]);
        }
        std::free(extef.cDx);
        std::free(extef.cDy);

        std::free(extef.bDx);
        std::free(extef.bDy);

        std::free (extef.intp);
    }

    void free (Vecteur& vect) { // inutile (mystère.)
        std::free (vect.toVector());
    }
}

#endif