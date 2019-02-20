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

void EIM::get_index(bool loud)
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

			if(loud) TLS.report(not_pol);
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

Rectangular::Rectangular(bool polarisation, double W, double H, double Ncore, double Nclad, double WL)
{
	// Primary Constructor

	set_params(polarisation, W, H, Ncore, Nclad, WL); 
}

void Rectangular::set_params(bool polarisation, double W, double H, double Ncore, double Nclad, double WL)
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
			pol = polarisation; 
			not_pol = !pol; 
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

void Rectangular::reduce_wg(bool loud)
{
	// reduce the 2D structure to two 1D structures

	try {
		if (params_defined) {
			//Step One: Eff Index along the vertical though the core

			TLS.set_params(height, lambda, ncore, nclad, nclad); // assign the parameters for the calculation

			TLS.neff_search(pol); // find the reduced effective index through the core

			if(loud) TLS.report(pol); 

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

// Wire Derived Class Definition
Wire::Wire()
{
	// Default Constructor
}

Wire::Wire(bool polarisation, double W, double H, double Ncore, double Nsub, double Nclad, double WL)
{
	set_params(polarisation, W, H, Ncore, Nsub, Nclad, WL);
}

void Wire::set_params(bool polarisation, double W, double H, double Ncore, double Nsub, double Nclad, double WL)
{
	// Assign values to the parameters of the rib waveguide class
	// R. Sheehan 20 - 2 - 2019

	try {
		bool c1 = W > 0.0 ? true : false;
		bool c2 = H > 0.0 ? true : false;
		bool c3 = Ncore > Nsub ? true : false;
		bool c4 = Nsub > Nclad ? true : false;
		bool c5 = Nclad > 0.0 ? true : false;
		bool c6 = WL > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5 && c6;

		if (c10) {
			pol = polarisation;
			not_pol = !pol;
			width = W; height = H;
			ncore = Ncore; nsub = Nsub; nclad = Nclad;
			lambda = WL;
			params_defined = true;
		}
		else {
			std::string reason;
			reason = "Error: void Wire::set_params(bool polarisation, double W, double H, double Ncore, double Nsub, double Nclad, double WL)\n";
			if (!c1) reason += "W: " + template_funcs::toString(W, 2) + " is not correct\n";
			if (!c2) reason += "H: " + template_funcs::toString(H, 2) + " is not correct\n";
			if (!c3) reason += "Ncore: " + template_funcs::toString(Ncore, 2) + " is not correct\n";
			if (!c4) reason += "Nsub: " + template_funcs::toString(Nsub, 2) + " is not correct\n";
			if (!c5) reason += "Nclad: " + template_funcs::toString(Nclad, 2) + " is not correct\n";
			if (!c6) reason += "WL: " + template_funcs::toString(WL, 2) + " is not correct\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void Wire::reduce_wg(bool loud)
{
	// reduce the 2D structure to two 1D structures

	try {
		if (params_defined) {
			//Step One: Eff Index along the vertical though the core

			TLS.set_params(height, lambda, ncore, nsub, nclad); // assign the parameters for the calculation

			TLS.neff_search(pol); // find the reduced effective index through the core

			if (loud) TLS.report(pol);

			if (TLS.get_nmodes(pol) > 0) {
				// Store the computed effective index values in the vectors eta_one, for cladding, and eta_two, for core
				for (int i = 0; i < TLS.get_nmodes(pol); i++) {
					eta_one.push_back(nclad);
					eta_two.push_back(TLS.get_neff(i, pol)); // store the reduced effective indices
				}
			}
			else {
				std::string reason;
				reason = "Error: void Wire::reduce_wg()\n";
				reason += "No modes computed in central vertical cross section\n";
				throw std::invalid_argument(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: void Wire::reduce_wg()\n";
			reason += "Device parameters not defined\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// Rib Derived Class Definition
Rib::Rib()
{
	// Default Constructor
}

Rib::Rib(bool polarisation, double W, double E, double T, double Ncore, double Nsub, double Nclad, double WL)
{
	set_params(polarisation, W, E, T, Ncore, Nsub, Nclad, WL); 
}

void Rib::set_params(bool polarisation, double W, double E, double T, double Ncore, double Nsub, double Nclad, double WL)
{
	// Assign values to the parameters of the rib waveguide class
	// R. Sheehan 20 - 2 - 2019

	try {
		bool c1 = W > 0.0 ? true : false;
		bool c2 = E > 0.0 ? true : false;
		bool c3 = T > 0.0 ? true : false;
		bool c4 = Ncore > Nsub ? true : false;
		bool c5 = Nsub > Nclad ? true : false;
		bool c6 = Nclad > 0.0 ? true : false;
		bool c7 = WL > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5 && c6 && c7;

		if (c10) {
			pol = polarisation;
			not_pol = !pol;
			width = W; etch_depth = E; slab_height = T; 
			ncore = Ncore; nsub = Nsub; nclad = Nclad;
			lambda = WL;
			params_defined = true;
		}
		else {
			std::string reason;
			reason = "Error: void Rib::set_params(bool polarisation, double W, double E, double T, double Ncore, double Nsub, double Nclad, double WL)\n";
			if (!c1) reason += "W: " + template_funcs::toString(W, 2) + " is not correct\n";
			if (!c2) reason += "E: " + template_funcs::toString(E, 2) + " is not correct\n";
			if (!c3) reason += "T: " + template_funcs::toString(T, 2) + " is not correct\n";
			if (!c4) reason += "Ncore: " + template_funcs::toString(Ncore, 2) + " is not correct\n";
			if (!c5) reason += "Nsub: " + template_funcs::toString(Nsub, 2) + " is not correct\n";
			if (!c6) reason += "Nclad: " + template_funcs::toString(Nclad, 2) + " is not correct\n";
			if (!c7) reason += "WL: " + template_funcs::toString(WL, 2) + " is not correct\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void Rib::reduce_wg(bool loud)
{
	
}

// Ridge Derived Class Definition
Shallow_Ridge::Shallow_Ridge()
{
	// Default Constructor
}

Shallow_Ridge::Shallow_Ridge(bool polarisation, double W, double E, double T, double D, double Ncore, double Nsub, double Nrib, double Nclad, double WL)
{
	set_params(polarisation, W, E, T, D, Ncore, Nsub, Nrib, Nclad, WL);
}

void Shallow_Ridge::set_params(bool polarisation, double W, double E, double T, double D, double Ncore, double Nsub, double Nrib, double Nclad, double WL)
{
	// Assign values to the parameters of the shallow ridge waveguide class
	// R. Sheehan 20 - 2 - 2019

	try {
		bool c1 = W > 0.0 ? true : false;
		bool c2 = E > 0.0 ? true : false;
		bool c3 = T > 0.0 ? true : false;
		bool c4 = D > 0.0 ? true : false;
		bool c5 = Ncore > Nrib ? true : false;
		bool c6 = Nsub > Nclad ? true : false;
		bool c7 = Nrib > Nclad ? true : false;
		bool c8 = Nclad > 0.0 ? true : false;
		bool c9 = WL > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8 && c9;

		if (c10) {
			pol = polarisation;
			not_pol = !pol;
			width = W; etch_depth = E; slab_height = T; core_height = D; 
			ncore = Ncore; nsub = Nsub; nrib = Nrib; nclad = Nclad;
			lambda = WL;
			params_defined = true;
		}
		else {
			std::string reason;
			reason = "Error: void Shallow_Ridge::set_params(bool polarisation, double W, double E, double T, double D, double Ncore, double Nsub, double Nrib, double Nclad, double WL)\n";
			if (!c1) reason += "W: " + template_funcs::toString(W, 2) + " is not correct\n";
			if (!c2) reason += "E: " + template_funcs::toString(E, 2) + " is not correct\n";
			if (!c3) reason += "T: " + template_funcs::toString(T, 2) + " is not correct\n";
			if (!c4) reason += "D: " + template_funcs::toString(D, 2) + " is not correct\n";
			if (!c5) reason += "Ncore: " + template_funcs::toString(Ncore, 2) + " is not correct\n";
			if (!c6) reason += "Nsub: " + template_funcs::toString(Nsub, 2) + " is not correct\n";
			if (!c7) reason += "Nrib: " + template_funcs::toString(Nrib, 2) + " is not correct\n";
			if (!c8) reason += "Nclad: " + template_funcs::toString(Nclad, 2) + " is not correct\n";
			if (!c9) reason += "WL: " + template_funcs::toString(WL, 2) + " is not correct\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void Shallow_Ridge::reduce_wg(bool loud)
{

}