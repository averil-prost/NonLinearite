#ifndef __INTERFACE__
#define __INTERFACE__

#include "inclusions.hpp"
#include "visu.hpp"

namespace interface {

    std::string entete_vtk (std::vector<Vecteur*> nodes, std::vector<int*> cells) {
        std::stringstream res;
        int i=0;
        // en-tête
        res << "# vtk DataFile Version 3.6" << std::endl;
        res << "PolyDATA" << std::endl;
        res << "ASCII" << std::endl;
        res << "DATASET UNSTRUCTURED_GRID" << std::endl;
        // points du maillage
        res << "POINTS\t" << nodes.size() << "\tfloat" << std::endl;
        for (Vecteur* vect : nodes) {
            for (i=0; i<3; ++i) {
                res << (*vect)(i) << (i<2?"\t":"\n");
            }
        }
        // cellules du maillage
        int ncells = cells.size();
        res << "CELLS\t" << ncells << "\t" << 4*ncells << std::endl; 
        for (int* cell : cells) {
            res << "3\t";
            for (i=1; i<4; ++i) {
                res << cell[i] << (i<3?"\t":"\n");
            }
        }
        // types de cellules
        res << "CELL_TYPES\t" << ncells << std::endl;
        for (i=0; i<ncells; ++i) {
            res << "5" << std::endl;
        }
        return res.str();
    }

    std::string vitesse_vtk (const Vecteur& Un, const Vecteur& dirich, int cardBaseV) {
        std::stringstream resstream;
        resstream << "POINT_DATA\t" << cardBaseV << std::endl;
        resstream << "VECTORS vitesse float" << std::endl;
        for (int i=0; i<cardBaseV; ++i) {
            resstream << Un(2*i)+dirich(2*i) << "\t" << Un(2*i+1)+dirich(2*i+1) << "\t0.0" << std::endl;
        }
        return resstream.str();
    }

    std::string pression_vtk (const Vecteur& Un, int cardBaseV, int cardBaseP) {
        std::stringstream resstream;
        resstream << "POINT_DATA\t" << cardBaseP << std::endl;
        resstream << "SCALARS pression float 1" << std::endl;
        resstream << "LOOKUP_TABLE default" << std::endl; 
        for (int i=0; i<cardBaseP; ++i) {
            resstream << Un(2*cardBaseV+i) << std::endl;
        }
        return resstream.str();
    }

    // filename sans arborescence ni extension
    void writeVTKAtStepN (const std::string& runname, int n, const std::string& entetev_vtk, const std::string& entetep_vtk, 
        const Vecteur& Un, const Vecteur& dirich, int cardBaseV, int cardBaseP) {
        std::ofstream vitesse_stream (VTK_FOLDER+runname+"_vitesse_"+std::to_string(n)+".vtk");
        vitesse_stream << entetev_vtk;
        vitesse_stream << vitesse_vtk (Un, dirich, cardBaseV); 
        vitesse_stream.close();

        std::ofstream pression_stream (VTK_FOLDER+runname+"_pression_"+std::to_string(n)+".vtk");
        pression_stream << entetep_vtk;
        pression_stream << pression_vtk (Un, cardBaseV, cardBaseP); 
        pression_stream.close();
    }

    // obsolete as paraview reads series just with name pattern
    void writeVTKtimeseries (const std::string& runname, int Nit, double deltaT) {
        std::ofstream vitesse_stream  (VTK_FOLDER+runname+"_vitesse_timeseries.series");
        std::ofstream pression_stream (VTK_FOLDER+runname+"_pression_timeseries.series");
        vitesse_stream  << "{\n\t\"file-series-version\" : \"1.0\",\n\t\"files\" : [\n";
        pression_stream << "{\n\t\"file-series-version\" : \"1.0\",\n\t\"files\" : [\n";
        for (int n=0; n<Nit; ++n) {
            vitesse_stream <<"\t\t{ \"name\" : \""<<VTK_FOLDER<<runname<<"_vitesse_" <<n<<".vtk\", \"time\" : "<<n*deltaT<<" }"<<(n==Nit-1?"":",")<<"\n";
            pression_stream<<"\t\t{ \"name\" : \""<<VTK_FOLDER<<runname<<"_pression_"<<n<<".vtk\", \"time\" : "<<n*deltaT<<" }"<<(n==Nit-1?"":",")<<"\n";
        }
        vitesse_stream  << "\t]\n}";
        pression_stream << "\t]\n}";
        vitesse_stream.close ();
        pression_stream.close();
    }

    std::string MatriceL_to_string (MatriceL& mat) {
        std::stringstream res;
        res << std::setprecision (14);
        int line=0, nelem=0, ielem=0, col=0; 
        res << mat.getDimX() << "\t" << mat.getDimY () << "\n";
        for (line = 0; line < mat.getDimX(); ++line) {
            nelem = mat.getNbTermL (line);
            res << nelem << (nelem==0 ? "\n" : "\t");
            for (ielem = 0; ielem < nelem; ++ielem) {
                col = mat.getIndLC (line, ielem);
                res << col << "\t" << mat (line, col) << (ielem==nelem-1?"\n":"\t");
            }
        }
        return res.str();
    }

    MatriceL MatriceL_from_string (const std::string& str) {
        std::stringstream stream;
        stream << str;
        int dimx=0, dimy=0, ielem=0, line=0, nelem=0, col=0;
        double value;
        stream >> dimx;
        stream >> dimy;
        MatriceL mat (dimx, dimy);
        for (line=0; line<dimx; ++line) {
            stream >> nelem;
            for (ielem=0; ielem<nelem; ++ielem) {
                stream >> col;
                stream >> value;
                mat.set (line, col, value);
            }
        }
        return mat;
    }

    void writeCube (const std::string which, const std::string mesh_name, int deg_ef_v, datastructure::Matrice3D& cube) {
        std::cout << "Writing " << which << std::endl;
        std::ofstream outstream ("data/"+which+"/"+which+"_"+mesh_name+"_deg"+std::to_string(deg_ef_v)+".txt");
        outstream << cube.cubeSize << "\n";
        int nnonz=0, nelem=0, ielem=0, col=0; 
        for (int k=0; k<cube.cubeSize; ++k) {
            nnonz = cube.nonzeros.at(k).size();
            outstream << nnonz << "\n"; // nombre de lignes j non nulles des matrices k 
            for (int line : cube.nonzeros.at(k)) {
                outstream << line << "\t"; // la ligne line est non nulle
                // copie de la ligne line de first
                nelem = cube.data.at(k).first.getNbTermL (line);
                outstream << nelem << "\t";
                for (ielem=0; ielem<nelem; ++ielem) {
                    col = cube.data.at(k).first.getIndLC (line, ielem);
                    outstream << col << "\t" << cube.data.at(k).first (line, col) << "\t";
                }
                // copie de la ligne line de second
                nelem = cube.data.at(k).second.getNbTermL (line);
                outstream << nelem << (nelem==0?"\n":"\t");
                for (ielem=0; ielem<nelem; ++ielem) {
                    col = cube.data.at(k).second.getIndLC (line, ielem);
                    outstream << col << "\t" << cube.data.at(k).second (line, col) << (ielem==nelem-1?"\n":"\t");
                }
            }
        } 
        outstream.close ();
        std::cout << which << " written." << std::endl;
    }

    datastructure::Matrice3D readCube (const std::string which, const std::string mesh_name, int deg_ef_v) {
        std::cout << "Loading " << which << std::endl;
        datastructure::Matrice3D cube;
        std::ifstream instream ("data/"+which+"/"+which+"_"+mesh_name+"_deg"+std::to_string(deg_ef_v)+".txt");
        int cubesize=0; instream >> cubesize; cube.cubeSize = cubesize;
        // plus simple de dupliquer que de convertir 14 fois des char* en stringstream en string...
        int k=0, ielem=0, line=0, nelem=0, col=0, nnonz=0, inonz=0;
        double value=0.0;
        for (k=0; k<cubesize; ++k) {
            MatriceL newmatx (cubesize), newmaty (cubesize);
            std::vector<int> nonz;
            instream >> nnonz;
            for (inonz=0; inonz<nnonz; ++inonz) {
                instream >> line;
                nonz.push_back (line);
                // copie de line pour first
                instream >> nelem;
                for (ielem=0; ielem<nelem; ++ielem) {
                    instream >> col;
                    instream >> value;
                    newmatx.set(line,col,value);
                }
                // copie de line pour second
                instream >> nelem;
                for (ielem=0; ielem<nelem; ++ielem) {
                    instream >> col;
                    instream >> value;
                    newmaty.set(line,col,value);
                }
            }
            cube.nonzeros.push_back (nonz);
            cube.data.push_back (std::make_pair<MatriceL, MatriceL> (std::move(newmatx), std::move(newmaty)));
        }

        instream.close ();
        std::cout << which << " loaded." << std::endl;
        return cube;
    }

    void writeMat (const std::string& which, const std::string mesh_name, int deg_ef_v, int deg_ef_p, MatriceL& mat) {
        std::cout << "Writing " << which << std::endl;
        std::ofstream outstream ("data/"+which+"/"+which+"_"+mesh_name+"_deg"+std::to_string(deg_ef_v)+(which=="B"?"-"+std::to_string(deg_ef_v):"")+".txt");
        outstream << interface::MatriceL_to_string (mat);
        outstream.close ();
        std::cout << which << " written." << std::endl;
    }

    MatriceL readMat (const std::string& which, const std::string mesh_name, int deg_ef_v, int deg_ef_p) {
        std::cout << "Loading " << which << std::endl;
        std::ifstream instream ("data/"+which+"/"+which+"_"+mesh_name+"_deg"+std::to_string(deg_ef_v)+(which=="B"?"-"+std::to_string(deg_ef_v):"")+".txt");
        std::stringstream stream;
        stream << instream.rdbuf(); // pfffffffff on se croirait en java
        MatriceL mat = interface::MatriceL_from_string (stream.str());
        instream.close ();
        std::cout << which << " loaded." << std::endl;
        return mat;
    }

    void writeVect (const std::string& which, const std::string& mesh_name, int deg, Vecteur& vect) {
        std::cout << "Writing " << which << std::endl;
        std::ofstream outstream ("data/"+which+"/"+which+"_"+mesh_name+"_deg"+std::to_string(deg)+".txt");
        int dim = vect.getDim();
        outstream << dim << "\n";
        for (int i=0; i<dim; ++i) {
            outstream << vect (i) << (i==dim-1?"\n":"\t");
        }
        outstream.close ();
        std::cout << which << " written." << std::endl;
    }

    Vecteur readVect (const std::string& which, const std::string& mesh_name, int deg) {
        std::cout << "Loading " << which << std::endl;
        std::ifstream instream ("data/"+which+"/"+which+"_"+mesh_name+"_deg"+std::to_string(deg)+".txt");
        int dim=0; 
        instream >> dim;
        double value=0.0;
        Vecteur res (dim);
        for (int i=0; i<dim; ++i) {
            instream >> value;
            res.set (i, value);
        }
        instream.close ();
        std::cout << which << " loaded." << std::endl;
        return res;
    }

    void exportSquareMesh (const std::string& mesh_name, double mesh_scale) {
        std::ofstream outstream (mesh_name+".geo");
        outstream << "Point(1) = {0, 0, 0, " << mesh_scale << "};" << std::endl;
        outstream << "Point(2) = {1, 0, 0, " << mesh_scale << "};" << std::endl;
        outstream << "Point(3) = {1, 1, 0, " << mesh_scale << "};" << std::endl;
        outstream << "Point(4) = {0, 1, 0, " << mesh_scale << "};" << std::endl;
        outstream << "Line(1) = {1, 2};" << std::endl;
        outstream << "Line(2) = {2, 3};" << std::endl;
        outstream << "Line(3) = {3, 4};" << std::endl;
        outstream << "Line(4) = {4, 1};" << std::endl;
        outstream << "Curve Loop(1) = {3, 4, 1, 2};" << std::endl;
        outstream << "Surface(1) = {1};" << std::endl;
        outstream.close ();

        system (("gmsh "+mesh_name+".geo -2 -o "+mesh_name+".msh").c_str());
    }

} // namespace interface

#endif