#ifndef ATTACH_H
#include "Attach.h"
#endif

// This is an implementation of the Effective Index Method for 2D waveguide cross sections
// It is valid for certain waveguide cross section geometries and is accurate to around 3 decimal places
// The method is based on the three and four layer slab waveguide solvers
// As a comparison consider the online solver http://www.computational-photonics.eu/eims.html
// R. Sheehan 13 - 4 - 2018

int main(int argc, char *argv[])
{
	//testing::slab_wg_mode_calc(); 

	//testing::eim_rect_wg(); 

	//testing::eim_wire_wg(); 

	//testing::eim_rib_wg(); 

	//testing::eim_ridge_wg(); 

	testing::eim_arb_wg(); 

	std::cout<<"Press enter to close console\n";
	std::cin.get(); 

	return 0; 
}