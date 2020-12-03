#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the methods declared in Eff_Indx_Method.h

// Definition of wg_dims class
wg_dims::wg_dims()
{
	// Default Constructor
	params_defined = false; 
	wg_type = 0; 
	width = height = etch_depth = slab_height = core_height = 0.0; 
}

wg_dims::wg_dims(wg_dims &obj)
{
	// Copy Constructor
	*this = obj; 
}

void wg_dims::set_dims(wg_dims &obj)
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
			width = W; height = H; etch_depth = slab_height = core_height = 0.0;
			wg_type = WIRE_WG; 
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
			width = W; etch_depth = E; slab_height = T; core_height = E + T; height = 0.0; 
			wg_type = RIB_WG; 
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
			width = W; etch_depth = E; slab_height = T; core_height = D; height = slab_height + etch_depth;
			wg_type = RIDGE_WG; 
			params_defined = true;
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

void ri_vals::set_ri(ri_vals &obj)
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
			ncore = Ncore; nsub = Nsub; nrib = Nrib; nclad = Nclad; lambda = WL; 
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
		bool c1 = dimensions.get_W() > 0.0 ? true : false; 
		bool c2 = eta_one.size() > 0 ? true : false; 
		bool c3 = eta_two.size() > 0 ? true : false;
		bool c4 = eta_two[0] > eta_one[0] ? true : false;
		bool c5 = ref_indx.get_WL() > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c5; 

		if (c10) {
			TLS.set_params(dimensions.get_W(), ref_indx.get_WL(), eta_two[0], eta_one[0], eta_one[0]);

			TLS.neff_search(not_pol); 

			if(loud) TLS.report(not_pol);
		}
		else {
			std::string reason;
			reason = "Error: void EIM::get_index()\n";
			if (!c1) reason += "width: " + template_funcs::toString(dimensions.get_W(), 2) + " is not correct\n";
			if (!c2) reason += "eta_one.size(): " + template_funcs::toString(eta_one.size(), 2) + " is not correct\n";
			if (!c3) reason += "eta_one.two(): " + template_funcs::toString(eta_two.size(), 2) + " is not correct\n";
			if (!c4) reason += "eta_two[0]: " + template_funcs::toString(eta_two[0], 2) + " <= eta_one[0]: " + template_funcs::toString(eta_one[0], 2) + " is not correct";
			if (!c5) reason += "lambda: " + template_funcs::toString(ref_indx.get_WL(), 2) + " is not correct\n";
			
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void EIM::set_params(bool polarisation, wg_dims &, ri_vals &)
{
	// Non-pure virtual function for the set_params function
	// Does nothing
}

double EIM::neff_value()
{
	// return the computed effective index value

	try {
		if(TLS.get_nmodes(not_pol) > 0){
			return TLS.get_neff(0, not_pol);
		}
		else {
			return 0.0; 
			std::string reason; 
			reason = "Error: double EIM::neff_value()\n"; 
			reason += "No modes computed during calculation\n"; 
			throw std::invalid_argument(reason); 
		}		
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

//double EIM::coupling_coefficient(double &separ)
//{
//	// estimate the coupling coefficient between this waveguide and a copy of itself
//	// theory taken from "Fundamentals of Optical Waveguides", Okamoto
//	// method uses the waveguide effective index to compute the approximation
//	// you could potentially vary the ridge width when doing this calculation, but adding another layer of objects to enable
//	// this would be very confusing
//	// calculation uses the info that was input into to enable calculation of 2D WG neff, so wavelength and material parameters are the same
//	// R. Sheehan 1 - 10 - 2018
//	// Added here 13 - 5 - 2019
//
//	try {
//		bool c1 = TLS.get_nmodes(not_pol) > 0 ? true : false; 
//		bool c2 = separ > 0 ? true : false;
//		bool c10 = c1 && c2; 
//		if (c10) {
//			return TLS.coupling_coefficient(separ, not_pol); 
//		}
//		else {
//			return 0.0;
//			std::string reason;
//			reason = "Error: double EIM::neff_value()\n";
//			if(!c1) reason += "No modes computed during calculation\n";
//			if (!c2) reason += "Waveguide separation must be positive\n"; 
//			throw std::invalid_argument(reason);
//		}
//	}
//	catch (std::invalid_argument &e) {
//		useful_funcs::exit_failure_output(e.what());
//		exit(EXIT_FAILURE);
//	}
//}

// Rectangular Derived Class Definition
Rectangular::Rectangular()
{
	// Default Constructor
}

Rectangular::Rectangular(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj)
{
	// Primary Constructor

	set_params(polarisation, dim_obj, ri_obj); 
}

void Rectangular::set_params()
{
	
}

//void Rectangular::set_params(bool polarisation, double W, double H, double Ncore, double Nclad, double WL)
//{
//	// Assign values to the parameters of the rectangular waveguide class
//	// R. Sheehan 20 - 2 - 2019
//
//	try {
//		bool c1 = W > 0.0 ? true : false; 
//		bool c2 = H > 0.0 ? true : false; 
//		bool c3 = Ncore > Nclad ? true : false; 
//		bool c4 = Nclad > 0.0 ? true : false; 
//		bool c5 = WL > 0.0 ? true : false; 
//		bool c10 = c1 && c2 && c3 && c4 && c5; 
//
//		if (c10) {
//			pol = polarisation; 
//			not_pol = !pol; 
//
//			width = W; height = H; 
//			
//			ncore = Ncore; nclad = Nclad; 
//			lambda = WL; 
//
//			dimensions.set_rect_wire(W, H); 
//
//			ref_indx.set_rect(Ncore, Nclad, WL); 
//
//			eta_one.clear(); 
//			eta_two.clear(); 
//
//			params_defined = true; 
//		}
//		else {
//			std::string reason; 
//			reason = "Error: void Rectangular::set_params(double W, double H, double Ncore, double Nclad, double WL)\n"; 
//			if (!c1) reason += "W: " + template_funcs::toString(W, 2) + " is not correct\n";
//			if (!c2) reason += "H: " + template_funcs::toString(H, 2) + " is not correct\n";
//			if (!c3) reason += "Ncore: " + template_funcs::toString(Ncore, 2) + " is not correct\n";
//			if (!c4) reason += "Nclad: " + template_funcs::toString(Nclad, 2) + " is not correct\n";
//			if (!c5) reason += "WL: " + template_funcs::toString(WL, 2) + " is not correct\n";
//			
//			throw std::invalid_argument(reason); 
//		}
//	}
//	catch (std::invalid_argument &e) {
//		useful_funcs::exit_failure_output(e.what());
//		exit(EXIT_FAILURE);
//	}
//}

void Rectangular::set_params(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj)
{
	// Non-pure virtual function for the set_params function
	// Assign values to the parameters of the rectangular waveguide class
	// R. Sheehan 28 - 2 - 2019

	try {
		bool c1 = dim_obj.get_W() > 0.0 ? true : false;
		bool c2 = dim_obj.get_H() > 0.0 ? true : false;
		bool c3 = ri_obj.get_Ncore() > ri_obj.get_Nclad() ? true : false;
		bool c4 = ri_obj.get_Nclad() > 0.0 ? true : false;
		bool c5 = ri_obj.get_WL() > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5;

		if (c10) {
			pol = polarisation;
			not_pol = !pol;

			dimensions.set_dims(dim_obj);

			ref_indx.set_ri(ri_obj);

			eta_one.clear();
			eta_two.clear();

			params_defined = true;
		}
		else {
			std::string reason;
			reason = "Error: void Rectangular::set_params(bool, wg_dims &, ri_vals &)\n";
			if (!c1) reason += "W: " + template_funcs::toString(dim_obj.get_W(), 2) + " is not correct\n";
			if (!c2) reason += "H: " + template_funcs::toString(dim_obj.get_H(), 2) + " is not correct\n";
			if (!c3) reason += "Ncore: " + template_funcs::toString(ri_obj.get_Ncore(), 2) + " is not correct\n";
			if (!c4) reason += "Nclad: " + template_funcs::toString(ri_obj.get_Nclad(), 2) + " is not correct\n";
			if (!c5) reason += "WL: " + template_funcs::toString(ri_obj.get_WL(), 2) + " is not correct\n";

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

			TLS.set_params(dimensions.get_H(), ref_indx.get_WL(), ref_indx.get_Ncore(), ref_indx.get_Nclad(), ref_indx.get_Nclad()); // assign the parameters for the calculation

			TLS.neff_search(pol); // find the reduced effective index through the core

			//if(loud) TLS.report(pol); 

			if (TLS.get_nmodes(pol) > 0) {
				// Store the computed effective index values in the vectors eta_one, for cladding, and eta_two, for core
				for (int i = 0; i < TLS.get_nmodes(pol); i++) {
					eta_one.push_back(ref_indx.get_Nclad());
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

Wire::Wire(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj)
{
	set_params(polarisation, dim_obj, ri_obj);
}

void Wire::set_params()
{

}

//void Wire::set_params(bool polarisation, double W, double H, double Ncore, double Nsub, double Nclad, double WL)
//{
//	// Assign values to the parameters of the rib waveguide class
//	// R. Sheehan 20 - 2 - 2019
//
//	try {
//		bool c1 = W > 0.0 ? true : false;
//		bool c2 = H > 0.0 ? true : false;
//		bool c3 = Ncore > Nsub ? true : false;
//		bool c4 = Nsub > Nclad ? true : false;
//		bool c5 = Nclad > 0.0 ? true : false;
//		bool c6 = WL > 0.0 ? true : false;
//		bool c10 = c1 && c2 && c3 && c4 && c5 && c6;
//
//		if (c10) {
//			pol = polarisation;
//			not_pol = !pol;
//			width = W; height = H;
//			ncore = Ncore; nsub = Nsub; nclad = Nclad;
//			lambda = WL;
//
//			dimensions.set_rect_wire(W, H); 
//
//			ref_indx.set_rib_wire(Ncore, Nsub, Nclad, WL); 
//
//			eta_one.clear();
//			eta_two.clear();
//			params_defined = true;
//		}
//		else {
//			std::string reason;
//			reason = "Error: void Wire::set_params(bool polarisation, double W, double H, double Ncore, double Nsub, double Nclad, double WL)\n";
//			if (!c1) reason += "W: " + template_funcs::toString(W, 2) + " is not correct\n";
//			if (!c2) reason += "H: " + template_funcs::toString(H, 2) + " is not correct\n";
//			if (!c3) reason += "Ncore: " + template_funcs::toString(Ncore, 2) + " is not correct\n";
//			if (!c4) reason += "Nsub: " + template_funcs::toString(Nsub, 2) + " is not correct\n";
//			if (!c5) reason += "Nclad: " + template_funcs::toString(Nclad, 2) + " is not correct\n";
//			if (!c6) reason += "WL: " + template_funcs::toString(WL, 2) + " is not correct\n";
//
//			throw std::invalid_argument(reason);
//		}
//	}
//	catch (std::invalid_argument &e) {
//		useful_funcs::exit_failure_output(e.what());
//		exit(EXIT_FAILURE);
//	}
//}

void Wire::set_params(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj)
{
	// Non-pure virtual function for the set_params function

	// Assign values to the parameters of the rib waveguide class
	// R. Sheehan 28 - 2 - 2019

	try {
		bool c1 = dim_obj.get_W() > 0.0 ? true : false;
		bool c2 = dim_obj.get_H() > 0.0 ? true : false;
		bool c3 = ri_obj.get_Ncore() > ri_obj.get_Nsub() ? true : false;
		bool c4 = ri_obj.get_Nsub() > ri_obj.get_Nclad() ? true : false;
		bool c5 = ri_obj.get_Nclad() > 0.0 ? true : false;
		bool c6 = ri_obj.get_WL() > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5 && c6;

		if (c10) {
			pol = polarisation;

			not_pol = !pol;
			
			dimensions.set_dims(dim_obj); 

			ref_indx.set_ri(ri_obj); 

			eta_one.clear();
			eta_two.clear();
			params_defined = true;
		}
		else {
			std::string reason;
			reason = "Error: void Wire::set_params(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj)\n";
			if (!c1) reason += "W: " + template_funcs::toString(dim_obj.get_W(), 2) + " is not correct\n";
			if (!c2) reason += "H: " + template_funcs::toString(dim_obj.get_H(), 2) + " is not correct\n";
			if (!c3) reason += "Ncore: " + template_funcs::toString(ri_obj.get_Ncore(), 2) + " is not correct\n";
			if (!c4) reason += "Nsub: " + template_funcs::toString(ri_obj.get_Nsub(), 2) + " is not correct\n";
			if (!c5) reason += "Nclad: " + template_funcs::toString(ri_obj.get_Nclad(), 2) + " is not correct\n";
			if (!c6) reason += "WL: " + template_funcs::toString(ri_obj.get_WL(), 2) + " is not correct\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void Wire::reduce_wg()
{
	// reduce the 2D structure to two 1D structures

	try {
		if (params_defined) {
			//Step One: Eff Index along the vertical though the core

			TLS.set_params(dimensions.get_H(), ref_indx.get_WL(), ref_indx.get_Ncore(), ref_indx.get_Nsub(), ref_indx.get_Nclad() ); // assign the parameters for the calculation

			TLS.neff_search(pol); // find the reduced effective index through the core

			//if (loud) TLS.report(pol);

			if (TLS.get_nmodes(pol) > 0) {
				// Store the computed effective index values in the vectors eta_one, for cladding, and eta_two, for core
				for (int i = 0; i < TLS.get_nmodes(pol); i++) {
					eta_one.push_back(ref_indx.get_Nclad());
					//eta_one.push_back(0.5*(nclad+nsub)); // Don't think this makes it more accurate
					eta_two.push_back(TLS.get_neff(i, pol)); // store the reduced effective indices
				}
			}
			else {
				// No modes through the vertical section means the device is basically a TLS anyway
				eta_one.push_back(ref_indx.get_Nclad());
				eta_two.push_back( (ref_indx.get_Ncore() + ref_indx.get_Nsub() + ref_indx.get_Nclad()) / 3.0);

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

Rib::Rib(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj)
{
	set_params(polarisation, dim_obj, ri_obj);
}

void Rib::set_params()
{

}

//void Rib::set_params(bool polarisation, double W, double E, double T, double Ncore, double Nsub, double Nclad, double WL)
//{
//	// Assign values to the parameters of the rib waveguide class
//	// R. Sheehan 20 - 2 - 2019
//
//	try {
//		bool c1 = W > 0.0 ? true : false;
//		bool c2 = E > 0.0 ? true : false;
//		bool c3 = T > 0.0 ? true : false;
//		bool c4 = Ncore > Nsub ? true : false;
//		bool c5 = Nsub > Nclad ? true : false;
//		bool c6 = Nclad > 0.0 ? true : false;
//		bool c7 = WL > 0.0 ? true : false;
//		bool c10 = c1 && c2 && c3 && c4 && c5 && c6 && c7;
//
//		if (c10) {
//			pol = polarisation;
//			not_pol = !pol;
//			width = W; etch_depth = E; slab_height = T; core_height = etch_depth + slab_height; 
//			ncore = Ncore; nsub = Nsub; nclad = Nclad;
//			lambda = WL;
//
//			dimensions.set_rib(W, E, T); 
//
//			ref_indx.set_rib_wire(Ncore, Nsub, Nclad, WL); 
//
//			eta_one.clear();
//			eta_two.clear();
//			params_defined = true;
//		}
//		else {
//			std::string reason;
//			reason = "Error: void Rib::set_params(bool polarisation, double W, double E, double T, double Ncore, double Nsub, double Nclad, double WL)\n";
//			if (!c1) reason += "W: " + template_funcs::toString(W, 2) + " is not correct\n";
//			if (!c2) reason += "E: " + template_funcs::toString(E, 2) + " is not correct\n";
//			if (!c3) reason += "T: " + template_funcs::toString(T, 2) + " is not correct\n";
//			if (!c4) reason += "Ncore: " + template_funcs::toString(Ncore, 2) + " is not correct\n";
//			if (!c5) reason += "Nsub: " + template_funcs::toString(Nsub, 2) + " is not correct\n";
//			if (!c6) reason += "Nclad: " + template_funcs::toString(Nclad, 2) + " is not correct\n";
//			if (!c7) reason += "WL: " + template_funcs::toString(WL, 2) + " is not correct\n";
//
//			throw std::invalid_argument(reason);
//		}
//	}
//	catch (std::invalid_argument &e) {
//		useful_funcs::exit_failure_output(e.what());
//		exit(EXIT_FAILURE);
//	}
//}

void Rib::set_params(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj)
{
	// Non-pure virtual function for the set_params function

	// Assign values to the parameters of the rib waveguide class
	// R. Sheehan 28 - 2 - 2019

	try {
		bool c1 = dim_obj.get_W() > 0.0 ? true : false;
		bool c2 = dim_obj.get_E() > 0.0 ? true : false;
		bool c3 = dim_obj.get_T() > 0.0 ? true : false;
		bool c4 = ri_obj.get_Ncore() > ri_obj.get_Nsub() ? true : false;
		bool c5 = ri_obj.get_Nsub() > ri_obj.get_Nclad() ? true : false;
		bool c6 = ri_obj.get_Nclad() > 0.0 ? true : false;
		bool c7 = ri_obj.get_WL() > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5 && c6 && c7;

		if (c10) {
			pol = polarisation;
			not_pol = !pol;
			
			dimensions.set_dims(dim_obj);

			ref_indx.set_ri(ri_obj);

			eta_one.clear();
			eta_two.clear();
			params_defined = true;
		}
		else {
			std::string reason;
			reason = "Error: void Rib::set_params(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj)\n";
			if (!c1) reason += "W: " + template_funcs::toString(dim_obj.get_W(), 2) + " is not correct\n";
			if (!c2) reason += "E: " + template_funcs::toString(dim_obj.get_E(), 2) + " is not correct\n";
			if (!c3) reason += "T: " + template_funcs::toString(dim_obj.get_T(), 2) + " is not correct\n";
			if (!c4) reason += "Ncore: " + template_funcs::toString(ri_obj.get_Ncore(), 2) + " is not correct\n";
			if (!c5) reason += "Nsub: " + template_funcs::toString(ri_obj.get_Nsub(), 2) + " is not correct\n";
			if (!c6) reason += "Nclad: " + template_funcs::toString(ri_obj.get_Nclad(), 2) + " is not correct\n";
			if (!c7) reason += "WL: " + template_funcs::toString(ri_obj.get_WL(), 2) + " is not correct\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void Rib::reduce_wg()
{
	// reduce the 2D structure to two 1D structures

	try {
		if (params_defined) {
			//Step One: Eff Index along the vertical though the slab

			TLS.set_params(dimensions.get_T(), ref_indx.get_WL(), ref_indx.get_Ncore(), ref_indx.get_Nsub(), ref_indx.get_Nclad() ); // assign the parameters for the calculation

			TLS.neff_search(pol); // find the reduced effective index through the core

			//if (loud) TLS.report(pol);

			if (TLS.get_nmodes(pol) > 0) {
				// Store the computed effective index values in the vectors eta_one, for slab, and eta_two, for core
				for (int i = 0; i < TLS.get_nmodes(pol); i++) {
					eta_one.push_back(TLS.get_neff(i, pol)); // store the reduced effective indices
				}
			}
			else {
				eta_one.push_back(ref_indx.get_Nclad());

				std::string reason;
				reason = "Error: void Rib::reduce_wg()\n";
				reason += "No modes computed in slab vertical cross section\n";
				throw std::runtime_error(reason);
			}

			//Step Two: Eff Index along the vertical though the core

			TLS.set_params(dimensions.get_D(), ref_indx.get_WL(), ref_indx.get_Ncore(), ref_indx.get_Nsub(), ref_indx.get_Nclad() ); // assign the parameters for the calculation

			TLS.neff_search(pol); // find the reduced effective index through the core

			//if (loud) TLS.report(pol);

			if (TLS.get_nmodes(pol) > 0) {
				// Store the computed effective index values in the vectors eta_one, for slab, and eta_two, for core
				for (int i = 0; i < TLS.get_nmodes(pol); i++) {
					eta_two.push_back(TLS.get_neff(i, pol)); // store the reduced effective indices
				}
			}
			else {
				eta_two.push_back((ref_indx.get_Ncore() + ref_indx.get_Nsub() + ref_indx.get_Nclad())/3.0);

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

Shallow_Ridge::Shallow_Ridge(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj)
{
	set_params(polarisation, dim_obj, ri_obj);
}

void Shallow_Ridge::set_params()
{

}

//void Shallow_Ridge::set_params(bool polarisation, double W, double E, double T, double D, double Ncore, double Nsub, double Nrib, double Nclad, double WL)
//{
//	// Assign values to the parameters of the shallow ridge waveguide class
//	// R. Sheehan 20 - 2 - 2019
//
//	try {
//		bool c1 = W > 0.0 ? true : false;
//		bool c2 = E > 0.0 ? true : false;
//		bool c3 = T > 0.0 ? true : false;
//		bool c4 = D > 0.0 ? true : false;
//		bool c5 = Ncore > Nrib ? true : false;
//		bool c6 = Nsub > Nclad ? true : false;
//		bool c7 = Nrib > Nclad ? true : false;
//		bool c8 = Nclad > 0.0 ? true : false;
//		bool c9 = WL > 0.0 ? true : false;
//		bool c10 = c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8 && c9;
//
//		if (c10) {
//			pol = polarisation;
//			not_pol = !pol;
//			width = W; etch_depth = E; slab_height = T; core_height = D; 
//			ncore = Ncore; nsub = Nsub; nrib = Nrib; nclad = Nclad;
//			lambda = WL;
//
//			dimensions.set_ridge(W, E, T, D); 
//
//			ref_indx.set_ridge(Ncore, Nsub, Nrib, Nclad, WL); 
//
//			eta_one.clear();
//			eta_two.clear();
//			params_defined = true;
//		}
//		else {
//			std::string reason;
//			reason = "Error: void Shallow_Ridge::set_params(bool polarisation, double W, double E, double T, double D, double Ncore, double Nsub, double Nrib, double Nclad, double WL)\n";
//			if (!c1) reason += "W: " + template_funcs::toString(W, 2) + " is not correct\n";
//			if (!c2) reason += "E: " + template_funcs::toString(E, 2) + " is not correct\n";
//			if (!c3) reason += "T: " + template_funcs::toString(T, 2) + " is not correct\n";
//			if (!c4) reason += "D: " + template_funcs::toString(D, 2) + " is not correct\n";
//			if (!c5) reason += "Ncore: " + template_funcs::toString(Ncore, 2) + " is not correct\n";
//			if (!c6) reason += "Nsub: " + template_funcs::toString(Nsub, 2) + " is not correct\n";
//			if (!c7) reason += "Nrib: " + template_funcs::toString(Nrib, 2) + " is not correct\n";
//			if (!c8) reason += "Nclad: " + template_funcs::toString(Nclad, 2) + " is not correct\n";
//			if (!c9) reason += "WL: " + template_funcs::toString(WL, 2) + " is not correct\n";
//
//			throw std::invalid_argument(reason);
//		}
//	}
//	catch (std::invalid_argument &e) {
//		useful_funcs::exit_failure_output(e.what());
//		exit(EXIT_FAILURE);
//	}
//}

void Shallow_Ridge::set_params(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj)
{
	// Non-pure virtual function for the set_params function

	// Assign values to the parameters of the shallow ridge waveguide class
	// R. Sheehan 28 - 2 - 2019

	try {
		bool c1 = dim_obj.get_W() > 0.0 ? true : false;
		bool c2 = dim_obj.get_E() > 0.0 ? true : false;
		bool c3 = dim_obj.get_T() > 0.0 ? true : false;
		bool c4 = dim_obj.get_D() > 0.0 ? true : false;
		bool c5 = ri_obj.get_Ncore() > ri_obj.get_Nrib() ? true : false;
		bool c6 = ri_obj.get_Nsub() > ri_obj.get_Nclad() ? true : false;
		bool c7 = ri_obj.get_Nrib() > ri_obj.get_Nclad() ? true : false;
		bool c8 = ri_obj.get_Nclad() > 0.0 ? true : false;
		bool c9 = ri_obj.get_WL() > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8 && c9;

		if (c10) {
			pol = polarisation;
			not_pol = !pol;			

			dimensions.set_dims(dim_obj);

			ref_indx.set_ri(ri_obj);

			eta_one.clear();
			eta_two.clear();
			params_defined = true;
		}
		else {
			std::string reason;
			reason = "Error: void Shallow_Ridge::set_params(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj)\n";
			if (!c1) reason += "W: " + template_funcs::toString(dim_obj.get_W(), 2) + " is not correct\n";
			if (!c2) reason += "E: " + template_funcs::toString(dim_obj.get_E(), 2) + " is not correct\n";
			if (!c3) reason += "T: " + template_funcs::toString(dim_obj.get_T(), 2) + " is not correct\n";
			if (!c4) reason += "D: " + template_funcs::toString(dim_obj.get_D(), 2) + " is not correct\n";
			if (!c5) reason += "Ncore: " + template_funcs::toString(ri_obj.get_Ncore(), 2) + " is not correct\n";
			if (!c6) reason += "Nsub: " + template_funcs::toString(ri_obj.get_Nsub(), 2) + " is not correct\n";
			if (!c7) reason += "Nrib: " + template_funcs::toString(ri_obj.get_Nrib(), 2) + " is not correct\n";
			if (!c8) reason += "Nclad: " + template_funcs::toString(ri_obj.get_Nclad(), 2) + " is not correct\n";
			if (!c9) reason += "WL: " + template_funcs::toString(ri_obj.get_WL(), 2) + " is not correct\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void Shallow_Ridge::reduce_wg()
{
	// reduce the 2D structure to two 1D structures

	try {
		if (params_defined) {
			//Step One: Eff Index along the vertical though the slab
			FLS.set_params(dimensions.get_D(), dimensions.get_T(), ref_indx.get_WL(), ref_indx.get_Ncore(), ref_indx.get_Nsub(), ref_indx.get_Nclad(), ref_indx.get_Nrib());
			
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

			FLS.set_params(dimensions.get_D(), dimensions.get_H(), ref_indx.get_WL(), ref_indx.get_Ncore(), ref_indx.get_Nsub(), ref_indx.get_Nclad(), ref_indx.get_Nrib());

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