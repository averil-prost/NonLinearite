// liste de toutes les inclusions "extérieures" 

#ifndef __EXTERN__
#define __EXTERN__

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <cstring>
#include <algorithm>
#include <set>

using namespace std; // grosse, grosse concession :'(

#include "../../../Lib/Vecteur.hpp"
#include "../../../Lib/MatriceL2.hpp"
#include "../../../Lib/element_fini_D1.hpp"
#include "../../../Lib/element_fini_D2.hpp"
#include "../../../Lib/element_fini_D3.hpp"
#include "../../../Lib/Meshes.hpp"
#include "../../../Lib/Bilinear.hpp"

#define MESH_FOLDER "data/gmsh/"
#define LINUX_SPECIFIC // comment if you are on anything but Linux
#define VTK_FOLDER "data/vtk/"

#endif