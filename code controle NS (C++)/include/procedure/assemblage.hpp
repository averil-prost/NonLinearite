///////////////////////////////////////////////////////////////////////////////////
//  partie assemblage du schéma
///////////////////////////////////////////////////////////////////////////////////

#ifndef __ASSEMBLAGE__
#define __ASSEMBLAGE__

#include "../inclusions.hpp"
#include "../datastructure.hpp"
#include "../builder.hpp"
#include "../interface.hpp"
#include "../visu.hpp"

namespace procedure {

    // forme trilinéaire 
    datastructure::Matrice3D assemblerCtilde (const datastructure::extension_ef_d2& extef,// const datastructure::bordermgr& bman,
        Meshes& mesh_v, const std::vector<int>& correspondance, int cardBaseV) {
        std::cout << "Getting ctilde" << std::endl;
        int phii=0,phij=0,phik=0;

        // initialisation de ctilde avec des matrices vides (tous les indices vont contenir qqch, pas besoin de sparse management (le vocable me manque))
        // std::cout << "Initialisation vide ctilde" << std::endl;
        datastructure::Matrice3D ctilde;
        ctilde.cubeSize = cardBaseV;
        for (int k=0; k<ctilde.cubeSize; ++k) {
            MatriceL matx (ctilde.cubeSize), maty (ctilde.cubeSize);
            ctilde.data.push_back (std::make_pair(std::move(matx),std::move(maty)));
        }

        // std::cout << "Itérations sur le maillage" << std::endl;
        std::vector<Vecteur*> nodes = mesh_v.getNodes();
        MatriceL JFmoins1 (2); // FOIS VALEUR ABSOLUE DU DÉTERMINANT
        double terme_dx=0.0, terme_dy=0.0;
        for (int* elem : mesh_v.getMesh2D()) { // pour chaque elem du maillage

            // DOCUMENTATION TIME les accès à ctilde
            //      notation théorique : ctilde_m (phi_i, phi_j, phi_k) avec m=1,2 indiquant par rapport à quelle variable on dérive
            //      notation du code : ctilde.at(k).<first/second> (j,i), first <=> m=1, second <=> m=2
            //   différents filtres pour savoir où placer quoi dans la construction : 
            //      phi_ijk : vont de 0 au nombre de fonctions de base de référence-1 (si P1, vont de 0 à 2), numérotation du treillis
            //      correspondance[^]+1 : dans le même intervalle, mais passe à la numérotation locale de gmsh 
            //      elem[^] : passage du local au global, donne le numéro de l'élément associé dans le maillage (converti en pk)
            //      bman.mesh_to_memory.at(^) : donne l'emplacement mémoire (même numérotation, mais dont on a viré les noeuds aux bords)
            // bonne chance au lecteur hihi

            // calcul des coefficients de la transformation : les 3 premiers noeuds donnent les sommets du triangle
            // inverse d'une matrice 2-2 : 1/det(matrice) * [[a22, -a12], [-a21, a11]]
            // en l'occurrence, inverse de F qui va de Kchapeau à K, donc annule le jacobien SAUF pour le signe
            JFmoins1.set (0,0,(*nodes[elem[3]])(1) - (*nodes[elem[1]])(1));
            JFmoins1.set (0,1,(*nodes[elem[1]])(0) - (*nodes[elem[3]])(0));
            JFmoins1.set (1,0,(*nodes[elem[1]])(1) - (*nodes[elem[2]])(1));
            JFmoins1.set (1,1,(*nodes[elem[2]])(0) - (*nodes[elem[1]])(0));
            if (JFmoins1(0,0)*JFmoins1(1,1)-JFmoins1(0,1)*JFmoins1(1,0)<0) JFmoins1 = JFmoins1 * -1.0; // signe de 1/det = signe de det

            // rangement : k est l'indice parcouru lors du parcours du vecteur
            // chaque matrice a j en ligne, i en colonne
            for (phii=0; phii<extef.cardBaseVref; ++phii) { // pour chaque fonction non nulle sur l'élément
                for (phij=0; phij<extef.cardBaseVref; ++phij) {
                    for (phik=0; phik<extef.cardBaseVref; ++phik) {
                        // std::cout << "\t\t(i,j,k) = (" << phii << "," << phij << "," << phik << ")\t" << std::endl;
                        terme_dx = 0.5*(extef.cDx[phii<=phik?phii:phik][phij][phii<=phik?phik-phii:phii-phik] - extef.cDx[phii<=phij?phii:phij][phik][phii<=phij?phij-phii:phii-phij]);
                        terme_dy = 0.5*(extef.cDy[phii<=phik?phii:phik][phij][phii<=phik?phik-phii:phii-phik] - extef.cDy[phii<=phij?phii:phij][phik][phii<=phij?phij-phii:phii-phij]);
                        // ok reader, don't freak out
                        // this is just ctilde.at(k).first(j,i) += JFmoins1(0,0)* terme_dx + JFmoins1(1,0) * terme_dy
                        // and the equivalent for ctilde.at(k).second
                        // see up for more details about how to access ctilde
                        ctilde.data.at(elem[correspondance[phik]+1]).first.set (elem[correspondance[phij]+1], elem[correspondance[phii]+1], 
                            ctilde.data.at(elem[correspondance[phik]+1]).first (elem[correspondance[phij]+1], elem[correspondance[phii]+1])
                            + JFmoins1(0,0) * terme_dx + JFmoins1(1,0) * terme_dy);
                        ctilde.data.at(elem[correspondance[phik]+1]).second.set (elem[correspondance[phij]+1], elem[correspondance[phii]+1], 
                            ctilde.data.at(elem[correspondance[phik]+1]).second (elem[correspondance[phij]+1], elem[correspondance[phii]+1])
                            + JFmoins1(0,1) * terme_dx + JFmoins1(1,1) * terme_dy);
                    }   
                }
            }
            // std::cout << std::endl;
        }

        // construction du vecteur des indices non vides
        // chaque matrice (pour chaque k) est non vide, mais à (k,j) fixé, il peut n'y avoir que des 0
        int line=0;
        for (auto kvalue : ctilde.data) { // pour chaque coordonnée k
            std::vector<int> nonzerosj;
            for (line=0; line<ctilde.cubeSize; ++line) { 
                if (kvalue.first.getNbTermL(line)>0 || kvalue.second.getNbTermL (line)>0) {
                    nonzerosj.push_back (line);
                }
            }
            ctilde.nonzeros.push_back (nonzerosj);
        }

        std::cout << "ctilde build." << std::endl;
        return ctilde;
    }

    // B dans le rapport associé 
    MatriceL assemblerMdivu1v2 (const datastructure::extension_ef_d2& extef, 
        Meshes& mesh_v, const std::vector<int>& correspondance_v, int cardBaseV,
        Meshes& mesh_p, const std::vector<int>& correspondance_p, int cardBaseP) {
        std::cout << "Getting mdivu1v2" << std::endl;

        MatriceL mdivu1v2 (2*cardBaseV, cardBaseP);
        int phii=0, psij=0;
        std::vector<Vecteur*> nodes_v = mesh_v.getNodes(), nodes_p = mesh_p.getNodes();
        MatriceL JFmoins1 (2); // FOIS VALEUR ABSOLUE DU DÉTERMINANT
        double terme_dx=0.0, terme_dy=0.0, value=0.0;
        
        int ielem=0, *elemp;
        std::vector<int*> pmesh2D = mesh_p.getMesh2D();
        for (int* elemv : mesh_v.getMesh2D()) { // pour chaque elem du maillage 
            // calcul des coefficients de la transformation : les 3 premiers noeuds donnent les sommets du triangle
            // inverse d'une matrice 2-2 : 1/det(matrice) * [[a22, -a12], [-a21, a11]]
            // en l'occurrence, inverse de F qui va de Kchapeau à K, donc annule le jacobien SAUF pour le signe
            JFmoins1.set (0,0,(*nodes_v[elemv[3]])(1) - (*nodes_v[elemv[1]])(1));
            JFmoins1.set (0,1,(*nodes_v[elemv[1]])(0) - (*nodes_v[elemv[3]])(0));
            JFmoins1.set (1,0,(*nodes_v[elemv[1]])(1) - (*nodes_v[elemv[2]])(1));
            JFmoins1.set (1,1,(*nodes_v[elemv[2]])(0) - (*nodes_v[elemv[1]])(0));
            if (JFmoins1(0,0)*JFmoins1(1,1)-JFmoins1(0,1)*JFmoins1(1,0)<0) JFmoins1 = JFmoins1 * -1.0; // signe de 1/det = signe de det

            elemp = pmesh2D.at(ielem);

            for (phii=0; phii<extef.cardBaseVref; ++phii) { // pour chaque fonction vitesse non nulle sur l'élément
                for (psij=0; psij<extef.cardBasePref; ++psij) { // pour chaque fonction pression non nulle sur l'élément
                    terme_dx = extef.bDx [phii][psij]; // ordre du treillis
                    terme_dy = extef.bDy [phii][psij];
                    mdivu1v2.set     (2*elemv[correspondance_v[phii]+1]  , elemp[correspondance_p[psij]+1], 
                            mdivu1v2 (2*elemv[correspondance_v[phii]+1]  , elemp[correspondance_p[psij]+1]) 
                            + JFmoins1(0,0)* terme_dx + JFmoins1(1,0) * terme_dy);
                    mdivu1v2.set     (2*elemv[correspondance_v[phii]+1]+1, elemp[correspondance_p[psij]+1], 
                            mdivu1v2 (2*elemv[correspondance_v[phii]+1]+1, elemp[correspondance_p[psij]+1]) 
                            + JFmoins1(0,1)* terme_dx + JFmoins1(1,1) * terme_dy);
                }
            }
            ++ielem;
        }
        std::cout << "mdivu1v2 build." << std::endl;
        return mdivu1v2;
    }

    // intégrales des fonctions de base de la pression
    Vecteur assemblerIntP (const datastructure::extension_ef_d2& extef,
        Meshes& mesh_p, const std::vector<int>& correspondance_p, int cardBaseP) {
        Vecteur intp (cardBaseP);
        int psij=0;
        double jacobien=0.0;
        std::vector<Vecteur*> nodes_p = mesh_p.getNodes();
        for (int* elem : mesh_p.getMesh2D()) { // pour chaque elem du maillage           
            jacobien = std::fabs(
                  ((*nodes_p[elem[2]])(0) - (*nodes_p[elem[1]])(0))*((*nodes_p[elem[3]])(1) - (*nodes_p[elem[1]])(1))
                - ((*nodes_p[elem[2]])(1) - (*nodes_p[elem[1]])(1))*((*nodes_p[elem[3]])(0) - (*nodes_p[elem[1]])(0))
            );
            for (psij=0; psij<extef.cardBasePref; ++psij) { // pour chaque fonction pression non nulle sur l'élément
                intp.set (elem[correspondance_p[psij]+1], intp(elem[correspondance_p[psij]+1]) + jacobien * extef.intp [psij]);
            }
        }
        return intp;
    }

    void tobeinvInsta (datastructure::insta& insta, const datastructure::bordermgr& bman,
            Vecteur intP, double nu, double deltaT) {
        std::cout << "Entrée dans tobeinvInsta" << std::endl;
        // principale utilité de cette fonction : masquer l'enfer de la copie
        
        // 1 terme en plus pour le multiplicateur de Lagrange de la moyenne nulle
        MatriceL tobeinverted (insta.cardBaseV*2 + insta.cardBaseP + 1);

        int ielem=0,col=0,nelem=0;
        int col_mass=0,ielem_mass=0,nelem_mass=0;
        int col_stif=0,ielem_stif=0,nelem_stif=0;
        double value=0.0;

        for (int line : bman.memory_to_mesh) { // copie des lignes associées aux points intérieurs
            // masse et rigidité (peut-être beaucoup de code pour rien...)
            nelem_mass = insta.uv.getNbTermL (line); 
            nelem_stif = insta.gradugradv.getNbTermL (line);
            ielem_mass = 0; ielem_stif = 0;
            while ((ielem_mass<nelem_mass) || (ielem_stif<nelem_stif)) {
                col_mass = insta.uv.getIndLC(line,ielem_mass); 
                col_stif = insta.gradugradv.getIndLC(line,ielem_stif); 
                if (col_mass == col_stif) {
                    value = insta.uv(line,col_mass)/deltaT + nu * insta.gradugradv(line,col_stif);
                    tobeinverted.set(2*line, 2*col_mass, value);
                    tobeinverted.set(2*line+1, 2*col_mass+1, value);
                    ielem_mass++;
                    ielem_stif++;
                }
                else if (col_mass < col_stif) {
                    value = insta.uv(line,col_mass)/deltaT;
                    tobeinverted.set(2*line, 2*col_mass, value);
                    tobeinverted.set(2*line+1, 2*col_mass+1, value);
                    ielem_mass++;
                }
                else {
                    value = nu * insta.gradugradv(line,col_stif);
                    tobeinverted.set(2*line, 2*col_stif, value);
                    tobeinverted.set(2*line+1, 2*col_stif+1, value);
                    ielem_stif++;
                }
            }
            // matrice des termes croisés 
            nelem = insta.mdivu1v2.getNbTermL (2*line);
            for (ielem=0; ielem<nelem; ++ielem) { // ligne pour la coord x
                col = insta.mdivu1v2.getIndLC (2*line, ielem);
                value = insta.mdivu1v2(2*line, col);
                tobeinverted.set (2*line, 2*insta.cardBaseV+col, value);
                tobeinverted.set (2*insta.cardBaseV+col, 2*line, value);
            }
            nelem = insta.mdivu1v2.getNbTermL (2*line+1);
            for (ielem=0; ielem<nelem; ++ielem) { // ligne pour la coord y
                col = insta.mdivu1v2.getIndLC (2*line+1, ielem);
                value = insta.mdivu1v2(2*line+1, col);
                tobeinverted.set (2*line+1, 2*insta.cardBaseV+col, value);
                tobeinverted.set (2*insta.cardBaseV+col, 2*line+1, value);
            }
        } // for line 

        // ajout de la condition de nullité de l'intégrale de P
        for (int line=0; line<insta.cardBaseP; ++line) {
            tobeinverted.set (2*insta.cardBaseV+insta.cardBaseP, 2*insta.cardBaseV+line, intP (line));
            tobeinverted.set (2*insta.cardBaseV+line, 2*insta.cardBaseV+insta.cardBaseP, intP (line));
        }

        insta.tobeinverted = tobeinverted;
        std::cout << "Sortie de tobeinvInsta" << std::endl;
    }

    MatriceL tobeinvSta (MatriceL& gradugradv, MatriceL& mdivu1v2, const datastructure::bordermgr& bman,
        Vecteur intP, double nu, int cardBaseV, int cardBaseP) {
        std::cout << "Entrée dans tobeinvSta" << std::endl;
        
        // 1 terme en plus pour le multiplicateur de Lagrange de la moyenne nulle
        MatriceL tobeinverted (cardBaseV*2 + cardBaseP + 1);

        int line=0,ielem=0,col=0,nelem=0;
        double value=0.0;
        for (line=0; line<cardBaseV; ++line) { 
            // rigidité
            nelem = gradugradv.getNbTermL (line);
            for (ielem=0; ielem<nelem; ++ielem) {
                col = gradugradv.getIndLC (line, ielem);
                value = nu * gradugradv(line,col);
                tobeinverted.set(2*line, 2*col, value);
                tobeinverted.set(2*line+1, 2*col+1, value);
            }
            // matrice des termes croisés 
            nelem = mdivu1v2.getNbTermL (2*line);
            for (ielem=0; ielem<nelem; ++ielem) { // ligne pour la coord x
                col = mdivu1v2.getIndLC (2*line, ielem);
                value = mdivu1v2(2*line, col);
                tobeinverted.set (2*line, 2*cardBaseV+col, value);
                tobeinverted.set (2*cardBaseV+col, 2*line, value);
            }
            nelem = mdivu1v2.getNbTermL (2*line+1);
            for (ielem=0; ielem<nelem; ++ielem) { // ligne pour la coord y
                col = mdivu1v2.getIndLC (2*line+1, ielem);
                value = mdivu1v2(2*line+1, col);
                tobeinverted.set (2*line+1, 2*cardBaseV+col, value);
                tobeinverted.set (2*cardBaseV+col, 2*line+1, value);
            }
        } // for line 

        // ajout de la condition de nullité de l'intégrale de P
        for (int line=0; line<cardBaseP; ++line) {
            tobeinverted.set (2*cardBaseV+cardBaseP, 2*cardBaseV+line, intP (line));
            tobeinverted.set (2*cardBaseV+line, 2*cardBaseV+cardBaseP, intP (line));
        }

        // std::cout << "Résultat de la copie : \n" << visu::tostring (tobeinverted);
        std::cout << "Sortie de tobeinvSta" << std::endl;
        return tobeinverted;
    }

    // calcule le second membre de stokes stationnaire pour t=0
    // namely       (f,v) - a(g0,v)     v \in V_h
    //              -b(g0,q)            q \in Y_h
    Vecteur initialScdMembre (MatriceL& uv, MatriceL& gradugradv, MatriceL& mdivu1v2, double nu, Vecteur& f0, Vecteur& dirich, 
        const datastructure::bordermgr& bman,  int cardBaseV, int cardBaseP) {
        
        std::cout << "Entrée dans initialScdMembre" << std::endl;
        Vecteur scdmembre (2*cardBaseV+cardBaseP+1);

        // M * f0 - A * g
        int nelem=0, ielem=0, col=0;
        double linevalue1=0.0, linevalue2=0.0, matvalue=0.0;
        for (int line : bman.memory_to_mesh) { // ne calculer que les noeuds intérieurs
            linevalue1 = 0.0; linevalue2=0.0;
            // M * f0
            nelem = uv.getNbTermL (line);
            for (ielem=0; ielem<nelem; ++ielem) {
                col = uv.getIndLC (line, ielem);
                matvalue = uv(line, col);
                linevalue1 += matvalue*f0(2*col);
                linevalue2 += matvalue*f0(2*col+1);
            }
            // - A * g
            nelem = gradugradv.getNbTermL (line);
            for (ielem=0; ielem<nelem; ++ielem) {
                col = gradugradv.getIndLC (line, ielem);
                matvalue = nu * gradugradv(line, col);
                linevalue1 -= matvalue*dirich(2*col);
                linevalue2 -= matvalue*dirich(2*col+1);
            }
            scdmembre.set (2*line, linevalue1);
            scdmembre.set (2*line+1, linevalue2);
        }

        // B^t * g
        // !!!!!! on dispose de B, on veut sa transposée => sparse dans le sens des colonnes
        // line parcourt les lignes de dirich
        double dirichlinevalue=0.0;
        for (int line=0; line<2*cardBaseV; ++line) {
            nelem = mdivu1v2.getNbTermL (line);
            dirichlinevalue = dirich (line);
            for (ielem=0; ielem<nelem; ++ielem) {
                col = mdivu1v2.getIndLC (line,ielem);
                scdmembre.set(2*cardBaseV+col, scdmembre(2*cardBaseV+col) - dirichlinevalue * mdivu1v2 (line, col));
            }
        }

        std::cout << "Sortie de initialScdMembre" << std::endl;
        return scdmembre;
    }

    // modifie les lignes associées au bord pour y mettre l'identité
    void pseudoElimination (MatriceL& matrix, const datastructure::bordermgr& bman, int cardBaseV) {
        std::cout << "Entrée dans pseudo-élimination" << std::endl;
        int nelem=0, ielem=0, col=0, jnode=0;
        double alpha = 0.0;
        bool stillInMat = false;
        
        // beaucoup de bruit pour rien (mais c'est fun) : ne pas perturber *du tout* le conditionnement
        // soit mat la partie de la matrice qui va rester, et nmat = [mat, 0; 0, alpha*Id] la nouvelle (à une permutation près)
        // plan : prendre alpha = racine d'une combinaison convexe des valeurs propres de mat^t * mat
        // comme ça, sp(nmat^t * nmat) = sp(mat^t * mat) U {alpha^2} et le conditionnement de nmat (en norme 2)
        // est inchangé. Donc on prend alpha = sqrt(trace(mat^t * mat)/taille(mat)), hihi
        for (int inode : bman.memory_to_mesh) { // pour chaque ligne de l'intérieur
            nelem = matrix.getNbTermL (2*inode); // partie x
            stillInMat = true; ielem=0;
            while (stillInMat && ielem < nelem) {
                col = matrix.getIndLC (2*inode, ielem);
                stillInMat = (col < 2*cardBaseV);
                if (stillInMat && bman.isInside.at(std::floor(col/2)))
                    alpha += matrix(2*inode, col) * matrix(col, 2*inode);
                ++ielem;
            }
            nelem = matrix.getNbTermL (2*inode+1); // partie y
            stillInMat = true; ielem=0;
            while (stillInMat && ielem < nelem) {
                col = matrix.getIndLC (2*inode+1, ielem);
                stillInMat = (col < 2*cardBaseV);
                if (stillInMat && bman.isInside.at(std::floor(col/2)))
                    alpha += matrix(2*inode+1, col) * matrix(col, 2*inode+1);
                ++ielem;
            }
        }
        alpha = std::sqrt(alpha / (2.0 * cardBaseV));

        int inode=0;
        for (int node=0; node < cardBaseV; ++node) {
            if (!bman.isInside.at(node)) {
                // première partie : les x
                nelem = matrix.getNbTermL (2*node);
                for (ielem=0; ielem<nelem; ++ielem) {
                    col = matrix.getIndLC (2*node, ielem);
                    matrix.set (2*node, col, 0.0);
                }
                // deuxième partie : les y
                nelem = matrix.getNbTermL (2*node+1);
                for (ielem=0; ielem<nelem; ++ielem) {
                    col = matrix.getIndLC (2*node+1, ielem);
                    matrix.set (2*node+1, col, 0.0);
                }
                // jusque-là, c'était propre, ça allait dans le sens des lignes
                // SUPPRESSION DES COLONNES => les fonctions associées ne sont pas
                // dans la base de l'espace test
                for (jnode=0; jnode<matrix.getDimY(); ++jnode) {
                    matrix.set (jnode, 2*node, 0.0);
                    matrix.set (jnode, 2*node+1, 0.0);
                }
                // ajout de l'identité x qqch
                matrix.set (2*node, 2*node, alpha);
                matrix.set (2*node+1, 2*node+1, alpha);
            }
        }
        std::cout << "Sortie de pseudo-élimination" << std::endl;
    }

    datastructure::schemas assemblerschemas (const std::string& insta_name, const std::string& mesh_name, 
        element_fini_D2& ef_v, element_fini_D2& ef_p, bool load_matrices, bool save_matrices,
        Meshes& mesh_v, Meshes& mesh_p, const datastructure::bordermgr& bman,
        double nu, double deltaT, int cardBaseV, int cardBaseP) {
        
        datastructure::schemas schemas;
        schemas.insta.name = insta_name;
        schemas.insta.cardBaseV = cardBaseV;
        schemas.insta.cardBaseP = cardBaseP;
        int deg_ef_v = ef_v.getOrdre ();
        int deg_ef_p = ef_p.getOrdre ();

        if (load_matrices) { // on considère qu'elles sont déjà construites !!!

            schemas.insta.uv = interface::readMat ("uv", mesh_name, deg_ef_v, deg_ef_p);
            schemas.insta.gradugradv = interface::readMat ("gradugradv", mesh_name, deg_ef_v, deg_ef_p);
            schemas.insta.mdivu1v2 = interface::readMat ("mdivu1v2", mesh_name, deg_ef_v, deg_ef_p);
            Vecteur intP = interface::readVect ("intP", mesh_name, deg_ef_p);

            procedure::tobeinvInsta (schemas.insta, bman, intP, nu, deltaT);
            if (insta_name != "stokes") {
                schemas.insta.ctilde = interface::readCube ("ctilde", mesh_name, deg_ef_v);
            }

            schemas.sta.tobeinverted = procedure::tobeinvSta (schemas.insta.gradugradv, schemas.insta.mdivu1v2, bman, intP, nu, cardBaseV, cardBaseP);

        } else { // construction (et enregistrement au passage) des matrices 

            // Le lecteur notera que c'est ici qu'est déterminé l'ordre des variables
            // on prend : (avec B = mdivu1v2)
            //  [           |      | 0 ] [v1_1]
            //  [  pas mal  |      | . ] [v1_2]
            //  [  de trucs |   B  | . ] [... ]     avec [vi_1, vi_2] la coordonnée de la vitesse
            //  [           |      | . ] [vN_1]     associée à la fonction phi_i
            //  [___________|______|_0_] [vN_2]     et le i dans l'ordre du maillage mesh_v
            //  [           |      |   ] [ p1 ]
            //  [    B^t    | rien | W ] [... ]
            //  [___________|______|___] [ pM ]
            //  [0   ...   0| W^t  | 0 ] [ λ  ]     multiplicateur de lagrange pour la pression

            // extension d'élément fini : contient les matrices élémentaires pour ctilde et mdivu1v2
            datastructure::extension_ef_d2 extef = builder::extension_ef_d2 (&ef_v, &ef_p);
            std::vector<int> correspondance_v = builder::get_corresp (*extef.ef_v, mesh_v); 
            std::vector<int> correspondance_p = builder::get_corresp (*extef.ef_p, mesh_p);

            schemas.insta.mdivu1v2 = procedure::assemblerMdivu1v2 (extef, mesh_v, correspondance_v, cardBaseV, mesh_p, correspondance_p, cardBaseP);
            if (insta_name != "stokes") {
                schemas.insta.ctilde = procedure::assemblerCtilde (extef, mesh_v, correspondance_v, cardBaseV);
            }

            // version pseudo-élimination
            std::cout << "Getting uv2D & gradugradv...";
            MatriceL uv (cardBaseV), gradugradv (cardBaseV);
            Bilinear myBilinear;
            myBilinear.uv2D (mesh_v.getNodes(), mesh_v.getMesh2D(), *extef.ef_v, &uv);
            myBilinear.gradugradv2D (mesh_v.getNodes(), mesh_v.getMesh2D(), *extef.ef_v, &gradugradv);
            schemas.insta.uv = uv;
            schemas.insta.gradugradv = gradugradv;
            std::cout << "Ok for uv2D & gradugradv." << std::endl;

            // porte bien son nom. Sert à imposer \int_\Omega P = 0
            Vecteur intP = procedure::assemblerIntP (extef, mesh_p, correspondance_p, cardBaseP);

            procedure::tobeinvInsta (schemas.insta, bman, intP, nu, deltaT);
            schemas.sta.tobeinverted = procedure::tobeinvSta (schemas.insta.gradugradv, schemas.insta.mdivu1v2, bman, intP, nu, cardBaseV, cardBaseP);

            // TANT QU'À FAIRE on enregistre tout ça
            if (save_matrices) {
                interface::writeVect ("intP", mesh_name, deg_ef_p, intP);
                interface::writeMat ("uv", mesh_name, deg_ef_v, deg_ef_p, schemas.insta.uv);
                interface::writeMat ("gradugradv", mesh_name, deg_ef_v, deg_ef_p, schemas.insta.gradugradv);
                interface::writeMat ("mdivu1v2", mesh_name, deg_ef_v, deg_ef_p, schemas.insta.mdivu1v2);
                if (insta_name != "stokes") {
                    interface::writeCube ("ctilde", mesh_name, deg_ef_v, schemas.insta.ctilde);
                }
            }

            // un peu de nettoyage
            builder::free (extef);
        }

        // dans tous les cas : préparation de la matrice à inverser pour la pseudo-élimination
        // on fait en sorte que les lignes concernées ne soient jamais changées, 
        // ni dans cette matrice, ni dans le second membre
        procedure::pseudoElimination (schemas.insta.tobeinverted, bman, cardBaseV);
        procedure::pseudoElimination (schemas.sta.tobeinverted, bman, cardBaseV);

        return schemas;
    }
}

#endif