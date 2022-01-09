///////////////////////////////////////////////////////////////////////////////////
//  partie test
///////////////////////////////////////////////////////////////////////////////////

#ifndef __TEST__
#define __TEST__

#include "../inclusions.hpp"
#include "../datastructure.hpp"
#include "../builder.hpp"
#include "../interface.hpp"
#include "../visu.hpp"

namespace procedure {

    namespace test {

        // différentes symétries de ctilde
        datastructure::Matrice3D testCube (int cubesize, const std::string& version) {
            datastructure::Matrice3D cube;
            cube.cubeSize = cubesize;
            int i=0,j=0,k=0;
            if (version=="dernier k full 1, reste 0") {
                for (i=0; i<cube.cubeSize-1; ++i) {
                    MatriceL newmatx (cube.cubeSize), newmaty(cube.cubeSize);
                    cube.data.push_back (std::make_pair<MatriceL,MatriceL>(std::move(newmatx),std::move(newmaty)));
                }
                MatriceL newmatx (cube.cubeSize), newmaty(cube.cubeSize);
                for (i=0; i<cube.cubeSize; ++i) {
                    for (j=0; j<cube.cubeSize; ++j) {
                        newmatx.set(i,j,1.0);
                        newmaty.set(i,j,1.0);
                    }
                }
                cube.data.push_back (std::make_pair<MatriceL,MatriceL>(std::move(newmatx),std::move(newmaty)));
            }
            else if (version=="flemme de nommer") {
                for (i=0; i<cube.cubeSize; ++i) {
                    MatriceL newmatx (cube.cubeSize), newmaty(cube.cubeSize);
                    for (j=0; j<cube.cubeSize; ++j) {
                        // for (k=0; k<cube.cubeSize; ++k) {
                            newmatx.set(j,k,(j==1?0.0:4.0));
                            newmaty.set(j,k,(j==1?0.0:1.0));
                        // }
                    }
                    cube.data.push_back (std::make_pair<MatriceL,MatriceL>(std::move(newmatx),std::move(newmaty)));
                }
            }
            else {
                std::cout << "Version " << version << " not found ??" << std::endl;
            }
            return cube;
        }

        // calcule la divergence "numérique"
        Vecteur Bproduct (MatriceL& mdivu1v2, Vecteur& XX, bool isVitesse) {
            int cardBaseV = mdivu1v2.getDimX()/2, cardBaseP = mdivu1v2.getDimY(), i=0,j=0;
            Vecteur res (isVitesse ? cardBaseP : 2*cardBaseV);

            if (isVitesse) {
                // checker la somme des colonnes
                // for (j=0; j<cardBaseP; ++j) {
                //     sum = 0.0;
                //     for (i=0; i<cardBaseV; ++i) {
                //     }
                //     std::cout << "Sum = " << sum << (std::fabs(sum)>1e-12 ? " NON NUUUUUUL":"") << std::endl;
                //     norm += sum*sum;
                // }
                // std::cout << "my norm : " << std::sqrt(norm) << std::endl;

                // faire le produit mdivu1v2^t * XX, attention à la position des termes
                double linevalue = 0.0;
                for (j=0; j<cardBaseP; ++j) {
                    linevalue = 0.0;
                    for (i=0; i<cardBaseV; ++i) {
                        linevalue += mdivu1v2(2*i,j) * XX(2*i) + mdivu1v2(2*i+1,j) * XX(2*i+1);
                    }
                    res.set (j, linevalue);
                }
            }
            return res;
        }

        // teste si P est dans l'espace des fonctions à intégrale nulle
        double checkPression (Vecteur& Un, MatriceL& tobeinverted, int cardBaseV, int cardBaseP) {
            double integrale = 0.0;
            int i=0;
            for (i=0; i<cardBaseP; ++i) {
                integrale += Un (2*cardBaseV+i) * tobeinverted (2*cardBaseV+cardBaseP, 2*cardBaseV+i);
            }
            // std::cout << "Intégrale de la pression : " << integrale << std::endl;
            return integrale;
        }

        void setCubeToOnes (datastructure::Matrice3D& ctilde, int cardBaseV) {
            for (int k=0; k<cardBaseV; ++k) {
                ctilde.nonzeros.at(k).clear ();
                for (int j=0; j<cardBaseV; ++j) {
                    for (int i=0; i<cardBaseV; ++i) {
                        ctilde.data.at(k).first.set(j,i,1.0);
                        ctilde.data.at(k).second.set(j,i,1.0);
                    }
                    ctilde.nonzeros.at(k).push_back (j);
                }
            }
        }

    } // namespace test
}

#endif