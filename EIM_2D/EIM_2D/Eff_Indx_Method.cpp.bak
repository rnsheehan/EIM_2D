#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the methods declared in Eff_Indx_Method.h

// Definition of wg_dims class
wg_dims::wg_dims()
{
	// Default Constructor
	params_defined = false; 
	width = height = etch_depth = slab_height = core_height = 0.0; 
}

wg_dims::wg_dims(wg_dims &obj)
{
	// Copy Constructor
	*this = obj; 
}

void wg_dims::set_rect_wire(double W, double H)
{
	// assign the values of the rect / wire waveguide dimensions

	try {
		bool c1 = W > 0.0 ? true : false;
		bool c2 = H > 0.0 ? true : false;
		bool c10 = c1 && c2;

		if (c10) {
			width = W; etch_depth = 0.0; slab_height = 0.0; core_height = H;
			params_defined = true; 
		}
		else {
			std::string reason;
			reason = "Error: void wg_dims::set_rect_wire(double W, double H)\n";
			if (!c1) reason += "W: " + template_funcs::toString(W, 2) + " is not correct\n";
			if (!c2) reason += "H: " + template_funcs::toString(H, 2) + " is not correct\n";
			
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void wg_dims::set_rib(double W, double E, double T)
{
	// assign the values of the rib waveguide dimensions

	try {
		bool c1 = W > 0.0 ? true : false;
		bool c2 = E > 0.0 ? true : false;
		bool c3 = T > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			width = W; etch_depth = E; slab_height = T; core_height = E + T;
			params_defined = true;
		}
		else {
			std::string reason;
			reason = "Error: void wg_dims::set_rib(double W, double E, double T)\n";
			if (!c1) reason += "W: " + template_funcs::toString(W, 2) + " is not correct\n";
			if (!c2) reason += "E: " + template_funcs::toString(E, 2) + " is not correct\n";
			if (!c3) reason += "T: " + template_funcs::toString(T, 2) + " is not correct\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void wg_dims::set_ridge(double W, double E, double T, double D)
{
	// assign the values of the ridge waveguide dimensions

	try {
		bool c1 = W > 0.0 ? true : false;
		bool c2 = E > 0.0 ? true : false;
		bool c3 = T > 0.0 ? true : false;
		bool c4 = D > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4; 

		if (c10) {
			width = W; etch_depth = E; slab_height = T; core_height = D;
		}
		else {
			std::string reason;
			reason = "Error: void wg_dims::set_ridge(double W, double E, double T, double D)\n";
			if (!c1) reason += "W: " + template_funcs::toString(W, 2) + " is not correct\n";
			if (!c2) reason += "E: " + template_funcs::toString(E, 2) + " is not correct\n";
			if (!c3) reason += "T: " + template_funcs::toString(T, 2) + " is not correct\n";
			if (!c4) reason += "D: " + template_funcs::toString(D, 2) + " is not correct\n";
			
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

/// Definition of ri_vals class
ri_vals::ri_vals()
{
	// Default Constructor
	params_defined = false;
	ncore = nsub = nrib = nclad = lambda = 0.0; 
}

ri_vals::ri_vals(ri_vals &obj)
{
	// Copy Constructor
	*this = obj; 
}

void ri_vals::set_rect(double Ncore, double Nclad, double WL)
{
	// Assign the RI values for a rect wg

	try {
		bool c1 = Ncore > Nclad ? true : false;
		bool c2 = Nclad > 0.0 ? true : false;
		bool c3 = WL > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			ncore = Ncore; nsub = 0.0; nrib = 0.0; nclad = Nclad; lambda = WL; 
			params_defined = true;
		}
		else {
			std::string reason;
			reason = "Error: void ri_vals::set_rect(double Ncore, double Nclad, double WL)\n";
			if (!c1) reason += "Ncore: " + template_funcs::toString(Ncore, 2) + " is not correct\n";
			if (!c2) reason += "Nclad: " + template_funcs::toString(Nclad, 2) + " is not correct\n";
			if (!c3) reason += "WL: " + template_funcs::toString(WL, 2) + " is not correct\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void ri_vals::set_rib_wire(double Ncore, double Nsub, double Nclad, double WL)
{
	// Assign the RI values for a rib / wire waveguide

	try {
		bool c1 = Ncore > Nsub ? true : false;
		bool c2 = Nsub > Nclad ? true : false;
		bool c3 = Nclad > 0.0 ? true : false;
		bool c4 = WL > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4;

		if (c10) {
			ncore = Ncore; nsub = Nsub; nrib = 0.0; nclad = Nclad; lambda = WL; 
			params_defined = true;
		}
		else {
			std::string reason;
			reason = "Error: void ri_vals::set_rib_wire(double Ncore, double Nsub, double Nclad)\n";
			if (!c1) reason += "Ncore: " + template_funcs::toString(Ncore, 2) + " is not correct\n";
			if (!c2) reason += "Nsub: " + template_funcs::toString(Nsub, 2) + " is not correct\n";
			if (!c3) reason += "Nclad: " + template_funcs::toString(Nclad, 2) + " is not correct\n";
			if (!c4) reason += "WL: " + template_funcs::toString(WL, 2) + " is not correct\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void ri_vals::set_ridge(double Ncore, double Nsub, double Nrib, double Nclad, double WL)
{
	// Assign the RI values for a ridge waveguide

	try {
		bool c1 = Ncore > Nrib ? true : false;
		bool c2 = Nsub > Nclad ? true : false;
		bool c3 = Nrib > Nclad ? true : false;
		bool c4 = Nclad > 0.0 ? true : false;
		bool c5 = WL > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5; 

		if (c10) {
			ncore = Ncore; nsub = Nsub; nrib = Nrib; nclad = Nclad; 
			params_defined = true;
		}
		else {
			std::string reason;
			reason = "Error: void ri_vals::set_ridge(double Ncore, double Nsub, double Nrib, double Nclad)\n";
			if (!c1) reason += "Ncore: " + template_funcs::toString(Ncore, 2) + " is not correct\n";
			if (!c2) reason += "Nsub: " + template_funcs::toString(Nsub, 2) + " is not correct\n";
			if (!c3) reason += "Nrib: " + template_funcs::toString(Nrib, 2) + " is not correct\n";
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
		//bool c4 = eta_two.size() == eta_one.size() ? true : false;
		bool c5 = lambda > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c5; 

		if (c10) {
			TLS.set_params(width, lambda, eta_two[0], eta_one[0], eta_one[0]); 

			TLS.neff_search(not_pol); 

			if(loud) TLS.report(not_pol);
		}
		else {
			std::string reason;
			reason = "Error: void EIM::get_index()\n";
			if (!c1) reason += "width: " + template_funcs::toString(width, 2) + " is not correct\n";
			if (!c2) reason += "eta_one.size(): " + template_funcs::toString(eta_one.size(), 2) + " is not correct\n";
			if (!c3) reason += "eta_one.two(): " + template_funcs::toString(eta_two.size(), 2) + " is not correct\n";
			if (!c5) reason += "lambda: " + template_funcs::toString(lambda, 2) + " is not correct\n";
			
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void EIM::set_params(bool, wg_dims &, ri_vals &)
{
	// Non-pure virtual function for the set_params function
	// Does nothing
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

void Rectangular::set_params()
{
	
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

			dimensions.set_rect_wire(W, H); 

			ref_indx.set_rect(Ncore, Nclad, WL); 

			eta_one.clear(); 
			eta_two.clear(); 

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

void Rectangular::set_params(bool, wg_dims &, ri_vals &)
{
	// Non-pure virtual function for the set_params function
}

void Rectangular::reduce_wg()
{
	// reduce the 2D structure to two 1D structures

	try {
		if (params_defined) {
			//Step One: Eff Index along the vertical though the core

			TLS.set_params(height, lambda, ncore, nclad, nclad); // assign the parameters for the calculation

			TLS.neff_search(pol); // find the reduced effective index through the core

			//if(loud) TLS.report(pol); 

			if (TLS.get_nmodes(pol) > 0) {
				// Store the computed effective index values in the vectors eta_one, for cladding, and eta_two, for core
				for (int i = 0; i < TLS.get_nmodes(pol); i++) {
					eta_one.push_back(nclad);
					eta_two.push_back(TLS.get_neff(i, pol)); // store the reduced effective indices
				}
			}
			else {
				// No modes through the vertical section means the device is basically a TLS anyway
				eta_one.push_back(nclad);
				eta_two.push_back(0.5*(ncore+nclad));

				std::string reason; 
				reason = "Error: void Rectangular::reduce_wg()\n";
				reason += "No modes computed in central vertical cross section\n"; 
				throw std::runtime_error(reason); 
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
	catch (std::runtime_error &e) {
		std::cerr << e.what(); 
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

void Wire::set_params()
{

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

			dimensions.set_rect_wire(W, H); 

			ref_indx.set_rib_wire(Ncore, Nsub, Nclad, WL); 

			eta_one.clear();
			eta_two.clear();
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

void Wire::set_params(bool, wg_dims &, ri_vals &)
{
	// Non-pure virtual function for the set_params function
}

void Wire::reduce_wg()
{
	// reduce the 2D structure to two 1D structures

	try {
		if (params_defined) {
			//Step One: Eff Index along the vertical though the core

			TLS.set_params(height, lambda, ncore, nsub, nclad); // assign the parameters for the calculation

			TLS.neff_search(pol); // find the reduced effective index through the core

			//if (loud) TLS.report(pol);

			if (TLS.get_nmodes(pol) > 0) {
				// Store the computed effective index values in the vectors eta_one, for cladding, and eta_two, for core
				for (int i = 0; i < TLS.get_nmodes(pol); i++) {
					eta_one.push_back(nclad);
					//eta_one.push_back(0.5*(nclad+nsub)); // Don't think this makes it more accurate
					eta_two.push_back(TLS.get_neff(i, pol)); // store the reduced effective indices
				}
			}
			else {
				// No modes through the vertical section means the device is basically a TLS anyway
				eta_one.push_back(nclad); 
				eta_two.push_back( (ncore + nsub + nclad) / 3.0);

				std::string reason;
				reason = "Error: void Wire::reduce_wg()\n";
				reason += "No modes computed in central vertical cross section\n";
				throw std::runtime_error(reason);
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
	catch (std::runtime_error &e) {
		std::cerr << e.what();
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

void Rib::set_params()
{

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
			width = W; etch_depth = E; slab_height = T; core_height = etch_depth + slab_height; 
			ncore = Ncore; nsub = Nsub; nclad = Nclad;
			lambda = WL;

			dimensions.set_rib(W, E, T); 

			ref_indx.set_rib_wire(Ncore, Nsub, Nclad, WL); 

			eta_one.clear();
			eta_two.clear();
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

void Rib::set_params(bool, wg_dims &, ri_vals &)
{
	// Non-pure virtual function for the set_params function
}

void Rib::reduce_wg()
{
	// reduce the 2D structure to two 1D structures

	try {
		if (params_defined) {
			//Step One: Eff Index along the vertical though the slab

			TLS.set_params(slab_height, lambda, ncore, nsub, nclad); // assign the parameters for the calculation

			TLS.neff_search(pol); // find the reduced effective index through the core

			//if (loud) TLS.report(pol);

			if (TLS.get_nmodes(pol) > 0) {
				// Store the computed effective index values in the vectors eta_one, for slab, and eta_two, for core
				for (int i = 0; i < TLS.get_nmodes(pol); i++) {
					eta_one.push_back(TLS.get_neff(i, pol)); // store the reduced effective indices
				}
			}
			else {
				eta_one.push_back(nclad); 

				std::string reason;
				reason = "Error: void Rib::reduce_wg()\n";
				reason += "No modes computed in slab vertical cross section\n";
				throw std::runtime_error(reason);
			}

			//Step Two: Eff Index along the vertical though the core

			TLS.set_params(core_height, lambda, ncore, nsub, nclad); // assign the parameters for the calculation

			TLS.neff_search(pol); // find the reduced effective index through the core

			//if (loud) TLS.report(pol);

			if (TLS.get_nmodes(pol) > 0) {
				// Store the computed effective index values in the vectors eta_one, for slab, and eta_two, for core
				for (int i = 0; i < TLS.get_nmodes(pol); i++) {
					eta_two.push_back(TLS.get_neff(i, pol)); // store the reduced effective indices
				}
			}
			else {
				eta_two.push_back((ncore + nsub + nclad)/3.0); 

				std::string reason;
				reason = "Error: void Rib::reduce_wg()\n";
				reason += "No modes computed in central vertical cross section\n";
				throw std::runtime_error(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: void Rib::reduce_wg()\n";
			reason += "Device parameters not defined\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what();
	}
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

void Shallow_Ridge::set_params()
{

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

			dimensions.set_ridge(W, E, T, D); 

			ref_indx.set_ridge(Ncore, Nsub, Nrib, Nclad, WL); 

			eta_one.clear();
			eta_two.clear();
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

void Shallow_Ridge::set_params(bool, wg_dims &, ri_vals &)
{
	// Non-pure virtual function for the set_params function
}

void Shallow_Ridge::reduce_wg()
{
	// reduce the 2D structure to two 1D structures

	try {
		if (params_defined) {
			//Step One: Eff Index along the vertical though the slab
			FLS.set_params(core_height, slab_height, lambda, ncore, nsub, nclad, nrib);
			
			FLS.neff_search(pol); 
			
			//if (loud) FLS.report(pol);	

			if (FLS.get_nmodes(pol) > 0) {
				// Store the computed effective index values in the vectors eta_one, for slab, and eta_two, for core
				for (int i = 0; i < FLS.get_nmodes(pol); i++) {
					eta_one.push_back(FLS.get_neff(i, pol)); // store the reduced effective indices
				}
			}
			else {
				eta_one.push_back(nclad);

				std::string reason;
				reason = "Error: void Shallow_Ridge::reduce_wg()\n";
				reason += "No modes computed in slab vertical cross section\n";
				throw std::runtime_error(reason);
			}

			//Step Two: Eff Index along the vertical though the core

			FLS.set_params(core_height, slab_height + etch_depth, lambda, ncore, nsub, nclad, nrib);

			FLS.neff_search(pol);

			//if (loud) FLS.report(pol);

			if (FLS.get_nmodes(pol) > 0) {
				// Store the computed effective index values in the vectors eta_one, for slab, and eta_two, for core
				for (int i = 0; i < FLS.get_nmodes(pol); i++) {
					eta_two.push_back(FLS.get_neff(i, pol)); // store the reduced effective indices
				}
			}
			else {
				eta_two.push_back((nclad+nrib+ncore+nsub)/4.0);

				std::string reason;
				reason = "Error: void Shallow_Ridge::reduce_wg()\n";
				reason += "No modes computed in central vertical cross section\n";
				throw std::runtime_error(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: void Shallow_Ridge::reduce_wg()\n";
			reason += "Device parameters not defined\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what();
	}
}