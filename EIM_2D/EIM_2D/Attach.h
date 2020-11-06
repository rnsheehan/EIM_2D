#ifndef ATTACH_H
#define ATTACH_H

// Be careful when trying to minimise the librairies that you include
// e.g. iostream includes cmath indirectly at least on MSVS
// https://stackoverflow.com/questions/29338108/is-cmath-or-math-h-really-needed-compiles-without-it
// If you are going to use a library function (or macro/template/whatever), it's up to you to include the correct header.
// Otherwise your program compiling correctly is simply an accident.

#include <cstdlib>
#include <iostream> // cout, cin, cerr
#include <iomanip> // setw, setprecision, time

#include <string>
#include <sstream>
#include <fstream>

// need these for directory manipulation
#include <direct.h>
#include <errno.h>

#include <cmath>
#include <vector>
#include <complex>

#include <algorithm>

// Constants
static const double p=(atan(1.0)); // pi / 4
static const double Two_PI=(8.0*p); // 2 pi
static const double PI=(4.0*p); // pi
static const double PI_2=(2.0*p); // pi / 2
static const double PI_3=((4.0/3.0)*p); // pi / 3
static const double PI_4=(p); // pi / 4
static const double PI_5=((4.0/5.0)*p); // pi / 5
static const double PI_6=((2.0/3.0)*p); // pi / 6 

static const bool TE=1;
static const bool TM=0;
static const bool Ex=1; // Ex propagation is equivalent to TE mode E^{x} polarisation => TM followed by TE
static const bool Ey=0; // Ey propagation is equivalent to TM mode E^{y} polarisation => TE followed by TM

// integer code to determine the particular waveguide object type
static const int RECT_WG = 5001; 
static const int WIRE_WG = 5002;
static const int RIB_WG = 5003;
static const int RIDGE_WG = 5004;

static const double EPS=(3.0e-12);

static const double SPEED_OF_LIGHT=(3.0e14); // Speed of light in microns per second
static const double EPSILON=(8.85e-18); // Permittivity of free space in Farads per micron
static const double MU=(12.566e-13); // Permeability of free space in Henrys per micron
static const double ETA=sqrt(MU/EPSILON); // Impedance of free space

static const int MAX_PATH_LENGTH = 250; // max. length for a directory in Windows OS

static const std::string empty_str = "";
static const std::string dottxt=".txt";

static const std::complex<double> zero(0.0, 0.0);
static const std::complex<double> eye(0.0, 1.0);
static const std::complex<double> one(1.0, 0.0);

#include "Templates.h"
#include "Useful.h"
#include "Vector_Utils.h"
#include "Slab_WG.h"
#include "Eff_Indx_Method.h"
#include "Test_Functions.h"

#endif