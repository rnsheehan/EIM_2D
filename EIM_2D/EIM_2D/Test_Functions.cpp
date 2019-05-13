#ifndef ATTACH_H
#include "Attach.h"
#endif

void testing::slab_wg_neff_calc()
{
	// Run a three layer slab waveguide neff calculation and check the results
	// R. Sheehan 2 - 5 - 2018

	double W = 1.5; 
	double WL = 1.55; 
	double Nc = 3.38; 
	double Ns = 3.17; 
	double Ncl = 1.0; 

	slab_tl_neff sl_obj; 

	sl_obj.set_params(W, WL, Nc, Ns, Ncl); 

	sl_obj.neff_search(TE); 

	std::cout << "Complete\n"; 
}

void testing::fl_slab_wg_neff_calc()
{
	// Compute the effective index in a four layer slab waveguide
	// R. Sheehan 5 - 10 - 2018

	double W, Wr, WL, Nc, Ns, Nr, Ncl; 

	// Case B: Field Oscillating in Core Only
	// For there to be a solution one has to have ncl < nm < nc, where nm = Max(nr,ns)
	//W = 0.5; Wr = 0.5; WL = 1.55; // n_{TM} = 3.274464
	//W = 0.5; Wr = 1.0; WL = 1.55; // n_{TM} = 3.27562
	//W = 0.6; Wr = 0.4; WL = 1.55; // n_{TM} = 3.29152
	W = 0.6; Wr = 0.9; WL = 1.55; // n_{TM} = 3.293
	Ns = 3.17; Nc = 3.38; Nr = 3.17; Ncl = 1.0;

	slab_fl_neff_B fl_obj;

	fl_obj.set_params(W, Wr, WL, Nc, Ns, Ncl, Nr); 

	fl_obj.neff_search(TE); 
}

void testing::fl_slab_wg_mode_calc()
{
	// Compute the mode profile in a four layer slab waveguide
	// R. Sheehan 20 - 2 - 2019

	double W, Wr, WL, Nc, Ns, Nr, Ncl;

	// Case B: Field Oscillating in Core Only
	// For there to be a solution one has to have ncl < nm < nc, where nm = Max(nr,ns)
	W = 0.5; Wr = 0.5; WL = 1.55; // n_{TM} = 3.274464
	//W = 0.5; Wr = 1.0; WL = 1.55; // n_{TM} = 3.27562
	//W = 0.6; Wr = 0.4; WL = 1.55; // n_{TM} = 3.29152
	//W = 0.6; Wr = 0.9; WL = 1.55; // n_{TM} = 3.293
	Ns = 3.17; Nc = 3.38; Nr = 3.17; Ncl = 1.0;
}

void testing::eim_rect_wg()
{
	// Approximate the effective index of a rectangular waveguide using EIM
	// R. Sheehan 20 - 2 - 2019

	// Confirmed from online simulator that method is accurate to within 1e-3
	// Seems to be more accurate when computing TE modes (starting from TM)
	// R. Sheehan 20 - 2 - 2019

	bool polarisation = TE; 
	double WL, W, H, nc, ncl; 

	WL = 1.55; 
	W = 2; H = 1; 
	nc = 3.38; ncl = 3.17; 

	wg_dims dim; 

	dim.set_rect_wire(W, H); 

	ri_vals ri; 

	ri.set_rect(nc, ncl, WL); 

	Rectangular wguide; 

	//wguide.set_params(polarisation, W, H, nc, ncl, WL);

	wguide.set_params(polarisation, dim, ri); 

	wguide.reduce_wg(); 

	wguide.get_index(true); 
}

void testing::eim_wire_wg()
{
	// Approximate the effective index of a wire waveguide using EIM
	// R. Sheehan 20 - 2 - 2019

	// Confirmed from online simulator that method is accurate to within 1e-3
	// Seems to be more accurate when computing TE modes (starting from TM)
	// When dealing with wire waveguides the TE mode (starting from TM) approximates 
	// the TM mode of the online simulator and vice versa
	// There must be a relationship between accuracy of EIM and index contrast
	// R. Sheehan 20 - 2 - 2019

	bool polarisation = TE;
	double WL, W, H, nc, ns, ncl;

	WL = 1.55;
	//W = 1; H = 0.5;
	//nc = 2.45; ns = 1.45; ncl = 1.0;
	// Si Wire has one TE mode with neff = 2.124
	// EIM starting with TE gives neff = 2.332 when ni is averaged over ns and ncl
	// EIM starting with TE gives neff = 2.294 when ni is ncl
	// EIM starting with TM gives neff = 1.625 when ni is averaged over ns and ncl
	// EIM starting with TM gives neff = 1.594 when ni is ncl
	W = 0.45; H = 0.22; 
	nc = 3.45; ns = 1.45; ncl = 1.0; 

	wg_dims dim;

	dim.set_rect_wire(W, H);

	ri_vals ri;

	ri.set_rib_wire(nc, ns, ncl, WL);

	Wire wguide;

	//wguide.set_params(polarisation, W, H, nc, ns, ncl, WL);

	wguide.set_params(polarisation, dim, ri); 

	wguide.reduce_wg();

	wguide.get_index(true);
}

void testing::eim_rib_wg()
{
	// Estimate the effective index of a rib waveguide
	// R. Sheehan 21 - 2 - 2019

	bool polarisation = TM; 

	double W, E, T, ncore, nsub, nclad, WL; 

	//W = 1.5; E = 0.3; T = 0.45; WL = 1.55; // neff = 3.265 for TM -> TE, neff = 3.269 for TE calc online
	//ncore = 3.38; nsub = 3.17; nclad = 1.0; // neff = 3.281 for TE -> TM, neff = 3.257 for TM calc online

	W = 2.0; E = 0.5; T = 0.5; WL = 1.55; // neff = 3.265 for TM -> TE, neff = 3.269 for TE calc online
	ncore = 3.38; nsub = 3.17; nclad = 1.0; // neff = 3.281 for TE -> TM, neff = 3.257 for TM calc online

	wg_dims dim;

	dim.set_rib(W, E, T);

	ri_vals ri;

	ri.set_rib_wire(ncore, nsub, nclad, WL);

	Rib wguide; 

	//wguide.set_params(polarisation, W, E, T, ncore, nsub, nclad, WL); 
	wguide.set_params(polarisation, dim, ri); 

	wguide.reduce_wg(); 

	wguide.get_index(true); 

	double separ = 0.25; 
	std::cout << separ << " , " << wguide.coupling_coefficient(separ) << "\n";
	separ = 0.5;
	std::cout << separ << " , " << wguide.coupling_coefficient(separ) << "\n";
	separ = 1.0; 
	std::cout << separ << " , " << wguide.coupling_coefficient(separ) << "\n";
	separ = 2.0;
	std::cout << separ << " , " << wguide.coupling_coefficient(separ) << "\n";
}

void testing::eim_ridge_wg()
{
	// Estimate the effective index of a ridge waveguide
	// R. Sheehan 21 - 2 - 2019

	bool polarisation = TM;

	double W, E, T, D, ncore, nsub, nrib, nclad, WL;

	W = 2.0; D = 0.6; E = 0.5; T = 0.4; WL = 1.55; // neff = 3.292 for TM -> TE, neff = 3.295 for TE calc online
	ncore = 3.38; nsub = nrib = 3.17; nclad = 1.0; 

	wg_dims dim;

	dim.set_ridge(W, E, T, D);

	ri_vals ri;

	ri.set_ridge(ncore, nsub, nrib, nclad, WL);

	Shallow_Ridge wguide;

	wguide.set_params(polarisation, dim, ri);

	wguide.reduce_wg();

	wguide.get_index(true);
}

void testing::eim_arb_wg()
{
	// Perform an EIM calculation using the EIM base class
	// R. Sheehan 28 - 2 - 2019

	bool polarisation = TM;

	double W, E, T, ncore, nsub, nclad, WL;

	W = 1.5; E = 0.3; T = 0.45; WL = 1.55; // neff = 3.265 for TM -> TE, neff = 3.269 for TE calc online
	ncore = 3.38; nsub = 3.17; nclad = 1.0; // neff = 3.281 for TE -> TM, neff = 3.257 for TM calc online

	//W = 2.0; E = 0.5; T = 0.5; WL = 1.55; // neff = 3.265 for TM -> TE, neff = 3.269 for TE calc online
	//ncore = 3.38; nsub = 3.17; nclad = 1.0; // neff = 3.281 for TE -> TM, neff = 3.257 for TM calc online

	wg_dims dim;

	dim.set_rib(W, E, T);

	ri_vals ri;

	ri.set_rib_wire(ncore, nsub, nclad, WL);

	Rib wguide; 

	EIM *compute = &wguide; 

	compute->set_params(polarisation, dim, ri); 

	compute->reduce_wg(); 

	compute->get_index(true); 
}