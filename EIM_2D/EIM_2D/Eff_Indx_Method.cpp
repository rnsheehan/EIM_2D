#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the methods declared in Eff_Indx_Method.h

// EIM Base Class Definition

EIM::EIM()
{
	// Default Constructor
	params_defined = false; 
	width = height = etch_depth = slab_height = core_height = 0.0; 
	ncore = nsub = nrib = nclad = lambda = 0.0; 
}

EIM::~EIM()
{
	// Deconstructor

	eta_one.clear(); 
	eta_two.clear(); 
	params_defined = false;
}

void EIM::get_index()
{
	// compute the effective index of the 2D structure
	// same calculation scheme will be used by all derived classes
	// R. Sheehan 31 - 5 - 2010
	try {
		bool c1 = width > 0.0 ? true : false; 
		bool c2 = eta_one.size() > 0 ? true : false; 
		bool c3 = eta_two.size() > 0 ? true : false;
		bool c4 = eta_two.size() == eta_one.size() ? true : false;
		bool c5 = lambda > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5; 

		if (c10) {
			TLS.set_params(width, lambda, eta_two[0], eta_one[0], eta_one[0]); 

			TLS.neff_search(not_pol); 
		}
		else {
			std::string reason;
			reason = "Error: void Rectangular::set_params(double W, double H, double Ncore, double Nclad, double WL)\n";
			if (!c1) reason += "width: " + template_funcs::toString(width, 2) + " is not correct\n";
			if (!c2 || c4) reason += "eta_one.size(): " + template_funcs::toString(eta_one.size(), 2) + " is not correct\n";
			if (!c3 || c4) reason += "eta_one.two(): " + template_funcs::toString(eta_two.size(), 2) + " is not correct\n";
			if (!c5) reason += "lambda: " + template_funcs::toString(lambda, 2) + " is not correct\n";
			
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// Rectangular Derived Class Definition
Rectangular::Rectangular()
{
	// Default Constructor
}

Rectangular::Rectangular(double W, double H, double Ncore, double Nclad, double WL)
{
	// Primary Constructor

	set_params(W, H, Ncore, Nclad, WL); 
}

void Rectangular::set_params(double W, double H, double Ncore, double Nclad, double WL)
{
	// Assign values to the parameters of the rectangular waveguide class
	// R. Sheehan 20 - 2 - 2019

	try {
		bool c1 = W > 0.0 ? true : false; 
		bool c2 = H > 0.0 ? true : false; 
		bool c3 = Ncore > Nclad ? true : false; 
		bool c4 = Nclad > 0.0 ? true : false; 
		bool c5 = WL > 0.0 ? true : false; 
		bool c10 = c1 && c2 && c3 && c4 && c5; 

		if (c10) {
			width = W; height = H; 
			ncore = Ncore; nclad = Nclad; 
			lambda = WL; 
			params_defined = true; 
		}
		else {
			std::string reason; 
			reason = "Error: void Rectangular::set_params(double W, double H, double Ncore, double Nclad, double WL)\n"; 
			if (!c1) reason += "W: " + template_funcs::toString(W, 2) + " is not correct\n";
			if (!c2) reason += "H: " + template_funcs::toString(H, 2) + " is not correct\n";
			if (!c3) reason += "Ncore: " + template_funcs::toString(Ncore, 2) + " is not correct\n";
			if (!c4) reason += "Nclad: " + template_funcs::toString(Nclad, 2) + " is not correct\n";
			if (!c5) reason += "WL: " + template_funcs::toString(WL, 2) + " is not correct\n";
			
			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void Rectangular::reduce_wg()
{
	// reduce the 2D structure to two 1D structures

	try {
		if (params_defined) {
			//Step One: Eff Index along the vertical though the core

			TLS.set_params(height, lambda, ncore, nclad, nclad); // assign the parameters for the calculation

			TLS.neff_search(pol); // find the reduced effective index through the core

			if (TLS.get_nmodes(pol) > 0) {
				// Store the computed effective index values in the vectors eta_one, for cladding, and eta_two, for core
				for (int i = 0; i < TLS.get_nmodes(pol); i++) {
					eta_one.push_back(nclad);
					eta_two.push_back(TLS.get_neff(i, pol)); // store the reduced effective indices
				}
			}
			else {
				std::string reason; 
				reason = "Error: void Rectangular::reduce_wg()\n";
				reason += "No modes computed in central vertical cross section\n"; 
				throw std::invalid_argument(reason); 
			}
		}
		else {
			std::string reason; 
			reason = "Error: void Rectangular::reduce_wg()\n"; 
			reason += "Device parameters not defined\n"; 

			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}