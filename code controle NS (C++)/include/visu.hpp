#ifndef __VISU__
#define __VISU__

#include "inclusions.hpp"

namespace visu {

    std::string tostring (MatriceL& mat) {
        std::stringstream res;
        res << std::scientific << std::setprecision(4);
        int i=0, j=0;
        for (i=0; i<mat.getDimX(); ++i) {
            for (j=0; j<mat.getDimY(); ++j) {
                // portable version on anything but Linux
                #ifndef LINUX_SPECIFIC
                res << mat(i,j) << "\t"; 
                #else
                // Linux-specific (and better) version
                if (std::fabs(mat(i,j)) <= 1e-10) {
                    res << "\033[2;37m";
                }
                res << std::setw (13) << mat(i,j);
                if (std::fabs(mat(i,j)) <= 1e-10) {
                    res << "\033[0m";
                }
                #endif
            }
            res << "\n";
        }
        return res.str();
    }

    std::string tostring (Vecteur& vect) {
        std::stringstream res;
        for (int i=0; i<vect.getDim(); ++i) {
            // portable version on anything but Linux
            #ifndef LINUX_SPECIFIC
            res << vect(i) << "\t"; 
            #else
            // Linux-specific (and better) version
            if (std::fabs(vect(i)) <= 1e-10) {
                res << "\033[2;37m";
            }
            res << std::setw (13) << vect(i);
            if (std::fabs(vect(i)) <= 1e-10) {
                res << "\033[0m";
            }
            #endif
        }
        return res.str()+"\n";
    }

    std::string tostring (Vecteur& vect, int cardBaseV, int cardBaseP) {
        std::stringstream res;
        res << std::scientific << std::setprecision(4);
        res << "Vitesse (" << cardBaseV << " éléments): " << std::endl;
        int i=0; 
        #ifndef LINUX_SPECIFIC
        for (i=0; i<cardBaseV; ++i) {
            res << "\t" << vect(2*i) << "\t" << vect(2*i+1) << std::endl;
        }
        res << "Pression (" << cardBaseP << " éléments): " << std::endl;
        for (i=0; i<cardBaseP; ++i) {
            res << "\t" << vect(2*cardBaseV+i) << std::endl;
        }
        #else
        for (i=0; i<cardBaseV; ++i) {
            if (std::fabs(vect(2*i)) <= 1e-10) res << "\033[2;37m";
            res << setw(12) << vect(2*i);
            if (std::fabs(vect(2*i)) <= 1e-10) res << "\033[0m";
            if (std::fabs(vect(2*i+1)) <= 1e-10) res << "\033[2;37m";
            res << setw(12) << vect(2*i+1) << std::endl;
            if (std::fabs(vect(2*i+1) <= 1e-10)) res << "\033[0m";
        }
        res << "Pression (" << cardBaseP << " éléments): " << std::endl;
        for (i=0; i<cardBaseP; ++i) {
            if (std::fabs(vect(2*cardBaseV+i)) <= 1e-10) res << "\033[2;37m";
            res << setw(12) << vect(2*cardBaseV+i) << std::endl;
            if (std::fabs(vect(2*cardBaseV+i)) <= 1e-10) res << "\033[0m";
        }
        #endif
        res << "Multiplicateur de Lagrange : " << vect (2*cardBaseV+cardBaseP) << std::endl;
        return res.str();
    }

    std::string tostring (double**tab, int dimx, int dimy) {
        std::stringstream res;
        res << std::scientific << std::setprecision(4);
        int i=0, j=0;
        for (i=0; i<dimx; ++i) {
            for (j=0; j<dimy; ++j) {
                // portable version on anything but Linux
                #ifndef LINUX_SPECIFIC
                res << tab[i][j] << "\t"; 
                #else
                // Linux-specific (and better) version
                if (std::fabs(tab[i][j]) <= 1e-10) {
                    res << "\033[2;37m";
                }
                res << std::setw (13) << tab[i][j];
                if (std::fabs(tab[i][j]) <= 1e-10) {
                    res << "\033[0m";
                }
                #endif
            }
            res << "\n";
        }
        return res.str();
    }

    std::string tostring (const datastructure::extension_ef_d2& extef) {
        std::stringstream res;
        res << "Extension d'un élément fini 2D d'ordre " << extef.ef_v->getOrdre() << std::endl;
        int i=0, j=0, k=0, m=1;
        // affichage de c
        for (m=1; m<=2; ++m) {
            res << "\nCube de la forme trilinéaire c de référence en d" << (m==1?"x":"y") << " : " << std::endl;
            for (k=0; k<extef.cardBaseVref; ++k) {
                res << "k = " << k << ", j ligne, i colonne -----------------------------------------" << std::endl;
                for (j=0; j<extef.cardBaseVref; ++j) {
                    for (i=0; i<extef.cardBaseVref; ++i) {
                    res << "\t";
                        res << std::setw (13) << (m==1?extef.cDx[i<=k?i:k][j][i<=k?k-i:i-k]:extef.cDy[i<=k?i:k][j][i<=k?k-i:i-k]) << "\t";
                    }
                    res << std::endl;
                }
            }
            res << std::endl;
        }
        res << std::endl;
        // affichage de ctilde
        for (m=1; m<=2; ++m) {
            res << "\nCube de la forme trilinéaire ctilde de référence en d" << (m==1?"x":"y") << " : " << std::endl;
            for (k=0; k<extef.cardBaseVref; ++k) {
                res << "k = " << k << ", j ligne, i colonne -----------------------------------------" << std::endl;
                for (j=0; j<extef.cardBaseVref; ++j) {
                    for (i=0; i<extef.cardBaseVref; ++i) {
                    res << "\t";
                        res << std::setw (13) << (m==1?0.5*(extef.cDx[i<=k?i:k][j][i<=k?k-i:i-k] - extef.cDx[i<=j?i:j][k][i<=j?j-i:i-j])
                                                      :0.5*(extef.cDy[i<=k?i:k][j][i<=k?k-i:i-k] - extef.cDy[i<=j?i:j][k][i<=j?j-i:i-j])) << "\t";
                    }
                    res << std::endl;
                }
            }
            res << std::endl;
        }
        res << std::endl;
        // affichage de b
        for (m=1; m<=2; ++m) {
            res << "\nMatrice B de référence en d" << (m==1?"x":"y") << " : " << std::endl;
            for (i=0; i<extef.cardBaseVref; ++i) {
                for (j=0; j<extef.cardBasePref; ++j) {
                    res << std::setw (13) << (m==1 ? extef.bDx[i][j] : extef.bDy[i][j]) << "\t";
                }
                res << std::endl;
            }
            res << std::endl;
        }
        // affichage des intégrales de la pression
        res << "\nIntégrales de la pression : ";
        for (j=0; j<extef.cardBasePref; ++j) {
            res << std::setw(13) << extef.intp [j] << "\t";
        }
        res << std::endl;
        return res.str();
    }

    std::string tostring (datastructure::Matrice3D& cube) {
        std::stringstream res;
        res << std::fixed;
        res << "Cube de dimension " << cube.cubeSize << std::endl;
        int k=0;
        for (k=0; k<cube.cubeSize; ++k) {
                res << "\nk = " << k << ", j ligne, i colonne -----------------------------------------" << std::endl;
                res << visu::tostring (cube.data.at(k).first) << std::endl;
                res << visu::tostring (cube.data.at(k).second);
                res << "Indices j supposément non nuls : ";
                for (int j : cube.nonzeros.at(k)) {
                    res << j << " ";
                }
                res << std::endl;
        }
        return res.str();
    }

    template <typename T> std::string tostring (const std::vector<T>& vect) {
        std::stringstream res;
        for (auto elem : vect) {
            res << std::setw(6) << elem << "\t";
        }
        res << std::endl;
        return res.str();
    }

    std::string tostring (poly_2D& poly) {
        std::stringstream res;
        int nbterms = poly.getNbTerm(), term=0;
        if (nbterms>0) res << poly.getMonome2D(0)->getCoeff();
        if (nbterms>1) {
            for (term=1; term<nbterms; ++term) {
                if (std::fabs(poly.getMonome2D(term)->getCoeff()) > 1e-8) {
                    res << (poly.getMonome2D(term)->getCoeff()>0?" + ":" - ") << \
                    (std::fabs(std::fabs(poly.getMonome2D(term)->getCoeff())-1.)<1e-8 ? "" : std::to_string(std::fabs(poly.getMonome2D(term)->getCoeff())));
                    res << (poly.getMonome2D(term)->getPowX()==0?"":"x" + (poly.getMonome2D(term)->getPowX()==1?"":"^"+std::to_string(poly.getMonome2D(term)->getPowX())) + " ");
                    res << (poly.getMonome2D(term)->getPowY()==0?"":"y" + (poly.getMonome2D(term)->getPowY()==1?"":"^"+std::to_string(poly.getMonome2D(term)->getPowY())));
                }
            }
        }
        res << std::endl;
        return res.str();
    }

    std::string tostring (Meshes& mesh) {
        std::stringstream res;
        // affichage des coordonnées des noeuds 
        int n=0;
        res << "Affichage du maillage" << std::endl;
        for (auto node : mesh.getNodes()) {
            res << "\tNoeud " << n++ << " : (" << (*node)(0) << ", " << (*node)(1) << ", " << (*node)(2) << ")" << std::endl; 
        }
        return res.str();
    }

    std::string tostring (const datastructure::bordermgr& bman) {
        std::stringstream res;
        int n=0;
        res << "Border manager : " << bman.nBorderNodes << " noeuds à la frontière" << std::endl;
        for (bool is : bman.isInside) {
            res << "\tNoeud " << n << (is ? " est à l'intérieur" : " est à la frontière");
            res << ", en mémoire sous " << bman.mesh_to_memory.at(n) << std::endl;
            n++;
        }
        res << "-------------------------------------------------------------------\n";
        n=0;
        for (int tomesh : bman.memory_to_mesh) {
            res << "\t" << n << "e emplacement mémoire conduit au noeud " << tomesh << std::endl;
            n++;
        }
        return res.str();
    }

    std::string stats (Vecteur& Un, int cardBaseV, int cardBaseP, double erreur_v, double erreur_p, double intP) {
        std::stringstream res;
        double moy_vitesse=0.0, moy_pression=0.0;
        int i=0;
        for (i=0; i<cardBaseV; ++i) {
            moy_vitesse += std::sqrt(Un(2*i)*Un(2*i) + Un(2*i+1)*Un(2*i+1));
        }
        moy_vitesse /= cardBaseV;
        for (i=0; i<cardBaseP; ++i) {
            moy_pression += Un (2*cardBaseV+i);
        }
        moy_pression /= cardBaseP;
        res << "Stats : magnitude moyenne vitesse = " << moy_vitesse;
        res << ", pression moyenne = " << moy_pression;
        res << ", lagrange = " << Un (2*cardBaseV+cardBaseP); 
        res << ", erreur vitesse = " << erreur_v;
        res << ", erreur pression = " << erreur_p;
        res << ", intégrale pression = " << intP;
        res << std::endl; 
        return res.str();
    }
}

#endif