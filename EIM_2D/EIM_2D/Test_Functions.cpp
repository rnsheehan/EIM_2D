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

void testing::slab_wg_mode_calc()
{
	// Run a three layer slab waveguide mode calculation and check the results
	// R. Sheehan 2 - 5 - 2018

	double W, WL, Nc, Ns; 

	//W = 1.0;  WL = 1.55; 
	//Nc = 3.38; Ns = 3.17; Ncl = 3.17;

	//W = 1.7;  WL = 6.158;
	//Nc = 3.54102488; Ns = 1.75275048;

	/*W = 2.0; WL = 1.55; 
	Nc = 3.27562; Ns = 3.27446;*/

	W = 2.5; WL = 1.55; 
	Nc = 3.293; Ns = 3.29152; 

	slab_tl_mode sl_obj;

	sl_obj.set_params(W, WL, Nc, Ns, Ns);

	sl_obj.compute_neff(TE); 

	sl_obj.compute_neff(TM);

	int N = 501; 
	double Lz = 3.0; 
	std::string stor = ""; 

	sl_obj.output_all_stats(stor);

	sl_obj.output_modes(TE, N, Lz, stor); 

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

	slab_fl_mode_B fl_obj;

	fl_obj.set_params(W, Wr, WL, Nc, Ns, Ncl, Nr); 

	fl_obj.compute_neff(TE); 
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

void testing::coupled_slab_wg_calc()
{
	// Compute the coupling coefficient for a pair of identical slab waveguides
	// R. Sheehan 1 - 10 - 2018

	double W, D, WL, Nc, Ns;

	// Code produces answer that matches with Okamoto example values
	// see Okamoto page 178, kappa = 0.39 mm^{-1}, L_{c} = 4 mm
	//W= 6; D = 12; WL = 1.5171; Nc = 1.46; Ns = 1.455;

	// Other examples

	// Si wire of height H = 0.22 um has neff^{TE} = 2.81 along vertical for nc = 3.45, ns = 1.45 and ncl = 1.0 at l = 1.55
	// Coupling coefficient between two Si wires of width W = 450 nm, D = 300 nm is kappa = 2.4 um^{-1}, L_{c} = 653 nm
	//W = 0.45; D = 0.3; WL = 1.55; Nc = 2.81; Ns = 1.0;

	double kappa;

	W = 4.0; WL = 6.09996577; Nc = 3.377661931; Ns = 3.148794556;

	coupled_slab_tl_neff sl_obj;

	D = 0.5; 

	for (int i = 0; i < 6; i++) {

		sl_obj.set_params(D, W, WL, Nc, Ns);

		kappa = sl_obj.compute_coupling_coeff(TM);

		std::cout << "Separation = " << D << " um\n";

		std::cout << "Coupling coefficient = " << kappa << " um^{-1}\n";

		std::cout << "Coupling length = " << (PI / (2.0*kappa)) << " um\n";

		std::cout << "\n\n";

		D += 1.0; 

	}
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

	Rectangular wguide; 

	wguide.set_params(polarisation, W, H, nc, ncl, WL); 

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
	// R. Sheehan 20 - 2 - 2019

	bool polarisation = TE;
	double WL, W, H, nc, ns, ncl;

	WL = 1.55;
	W = 1; H = 0.5;
	nc = 2.45; ns = 1.45; ncl = 1.0;

	Wire wguide;

	wguide.set_params(polarisation, W, H, nc, ns, ncl, WL);

	wguide.reduce_wg();

	wguide.get_index(true);
}