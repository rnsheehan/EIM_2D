#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the methods declared in Coupled_WG.h
// R. Sheehan 4 - 12 - 2020

coupWG::coupWG()
{
	// Default constructor
	WGA = NULL; 
	WGB = NULL; 
	params_defined = false;
	waveguides_reduced = false; 
}

coupWG::coupWG(EIM* wgobjA, EIM* wgobjB)
{
	// primary constructor
	set_params(wgobjA, wgobjB); 
}

coupWG::coupWG(coupWG& obj) 
{
	// copy constructor
	*this = obj; 
}

coupWG::~coupWG()
{
	// Destructor
	WGA = NULL;
	WGB = NULL;
	params_defined = false; 
	waveguides_reduced = false; 
}

void coupWG::set_params(EIM* wgobjA, EIM* wgobjB, bool loud)
{
	// method for assigning the parameter values
	// all parameters in the waveguide objects must be assigned correctly before proceeding with calculation
	// the way this calculating is being set up means that you can compute coupling between arbitrary waveguide objects
	// in some cases this may produce physically meaningless results, but the calculation is possible
	// R. Sheehan 4 -  12 - 2020

	try {
		bool c1 = wgobjA->defined(); 
		bool c2 = wgobjB->defined(); 
		bool c10 = c1 && c2; 

		if (c10) {
			// Parameters in each waveguide object have to be defined before proceeding with calculation
			WGA = wgobjA; 

			WGB = wgobjB; 

			params_defined = true;
		}
		else {
			std::string reason;
			reason = "Error: void coupWG::set_params(EIM& wgobjA, EIM& wgobjB, bool loud)\n";
			if (!c1) reason += "Parameters not defined in wgobjA\n";
			if (!c2) reason += "Parameters not defined in wgobjB\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void coupWG::reduce_wg()
{
	// Perform the EIM calculations to reduce the 2D WG to 1D slabs
	// R. Sheehan 4 - 12 - 2020

	try {
		if (params_defined) {
			
			WGA->reduce_wg(); // reduce WGA 2D -> 1D

			WGB->reduce_wg(); // reduce WGB 2D -> 1D

			// check if the waveguides entered actually allow bound modes
			// to be fair if either waveguide does not allow bound modes 
			// an error will already have been thrown
			if ( !(WGA->neff_value() > 0.0) || !(WGB->neff_value() > 0.0) ) {
				waveguides_reduced = false; 

				std::string reason;
				reason = "Error: void coupWG::reduce_wg()\n";
				reason += "One of the waveguides has no bound modes\n";

				throw std::runtime_error(reason);
			}
			else {
				waveguides_reduced = true; 
			}
		}
		else {
			std::string reason;
			reason = "Error: void coupWG::reduce_wg()\n";
			reason += "Parameters not defined in waveguide objects\n";			

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error& e) {
		std::cerr << e.what();
	}
}

void coupWG::coupling_coeffs(double pitch, bool loud)
{
	// compute the coupling coefficients for two reduced waveguides

	try {
		if (waveguides_reduced) {
			// both waveguides reduced 2D -> 1D, coupling calculation can proceed

			//WGA->get_dims().get_W(); 
			//WGB->get_dims().get_W(); 
			//WGA->get_RI().get_WL(); 
			//WGA->reduced_core()[0]; 
			//WGB->reduced_core()[0]; 
			//WGA->reduced_cladding()[0]; 

			// instantiate the coupled_slabs object
			coupling.set_params(WGA->get_dims().get_W(), WGB->get_dims().get_W(), WGA->get_RI().get_WL(), 
				WGA->reduced_core()[0], WGB->reduced_core()[0], WGA->reduced_cladding()[0]); 

			coupling.compute_coefficients(pitch, loud); 
		}
		else {
			std::string reason;
			reason = "Error: void coupWG::coupling_coeffs()\n";
			reason += "Coupling coefficient calculation not possible\n";

			throw std::runtime_error(reason);
		}
	}
	catch (std::runtime_error& e) {
		std::cerr << e.what();
	}
}