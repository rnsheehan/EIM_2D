#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definitions for the slab base class

slab::slab(void)
{
	// Default Constructor
	M = 0;

	k_sqr = nc_sqr = ncl_sqr = ns_sqr = nr_sqr = nm_sqr = k_sqr_nc_sqr = k_sqr_ns_sqr = k_sqr_ncl_sqr = k_sqr_nr_sqr = 0.0;
	k_sqr_nm_sqr = aa = bb = d = tt = dr = l = nc = ns = ncl = nr = nm = etacr = etacs = etarcl = 0.0;
	V = na = k = lower = upper = upper = w = efieldint = hfieldint = 0.0;
}

slab::~slab(void)
{
	// Deconstructor
	betaE.clear();

	betaH.clear(); 
}

int slab::nbeta(bool mode)
{
	// return the number of computed propagation constants for a given polarisation
	// it can be the case that the predicted value M is greater than the actual number of modes in a waveguide

	return static_cast<int>(mode ? betaE.size() : betaH.size() );
}

double slab::h(int i, bool t)
{
	//Wavenumber in Core

	try {

		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl_mode::h(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt(k_sqr_nc_sqr - template_funcs::DSQR(betaE[i]));
				}
				else {
					return sqrt(k_sqr_nc_sqr - template_funcs::DSQR(betaH[i]));
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab::h(int i, bool t)\nNo modes have been computed\n");
			return 0; 
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab::p(int i, bool t)
{
	//Wavenumber in Substrate

	try {
		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab::p(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt(template_funcs::DSQR(betaE[i]) - k_sqr_ns_sqr);
				}
				else {
					return sqrt(template_funcs::DSQR(betaH[i]) - k_sqr_ns_sqr);
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab::p(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab::q(int i, bool t)
{
	//Wavenumber in Cladding

	try {

		if (nbeta(t) > 0) {

			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab::q(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt(template_funcs::DSQR(betaE[i]) - k_sqr_ncl_sqr);
				}
				else {
					return sqrt(template_funcs::DSQR(betaH[i]) - k_sqr_ncl_sqr);
				}
			}

		}
		else {
			throw std::invalid_argument("Error: double slab::q(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab::rA(int i, bool t)
{
	//Wavenumber in Rib for Case A of FL slab

	try {

		if (nbeta(t) > 0) {

			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab::rA(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt(k_sqr_nr_sqr - template_funcs::DSQR(betaE[i]));
				}
				else {
					return sqrt(k_sqr_nr_sqr - template_funcs::DSQR(betaH[i]));
				}
			}

		}
		else {
			throw std::invalid_argument("Error: double slab::r(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab::rB(int i, bool t)
{
	//Wavenumber in Rib for Case B of FL slab

	try {

		if (nbeta(t) > 0) {

			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab::rB(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt(template_funcs::DSQR(betaE[i]) - k_sqr_nr_sqr);
				}
				else {
					return sqrt(template_funcs::DSQR(betaH[i]) - k_sqr_nr_sqr);
				}
			}

		}
		else {
			throw std::invalid_argument("Error: double slab::r(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error& e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab::_beta(int i, bool t)
{
	// return the i^{th} propagation constant for a given polarisation
	// by convention t == true => TE modes, t == false => TM modes

	try {
		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_wg::_beta(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return betaE[i];
				}
				else {
					return betaH[i];
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab::_beta(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// Definitions for the three layer slab derived class
slab_tl_neff::slab_tl_neff(void)
{
	// Default Constructor
	g = 0; 
}

slab_tl_neff::slab_tl_neff(double width, double lambda, double ncore, double nsub, double nclad)
{
	set_params(width, lambda, ncore, nsub, nclad);
}

void slab_tl_neff::set_params(double width, double lambda, double ncore, double nsub, double nclad)
{
	// assign values to the parameters for the three layer slab waveguide
	// throw exception if one of the inputs is not correct

	try {

		bool c1 = width > 0.0 ? true : false;
		bool c2 = lambda > 0.0 ? true : false;
		bool c3 = nclad >= 1.0 ? true : false;
		bool c4 = nsub >= 1.0 ? true : false;
		bool c5 = ncore > std::max(nsub, nclad) ? true : false;

		if (c1 && c2 && c3 && c4 && c5) {

			d = width;

			tt = 0.5 * d; 

			l = lambda;

			nc = ncore;
			nc_sqr = template_funcs::DSQR(nc);

			ns = std::max(nsub, nclad);
			ns_sqr = template_funcs::DSQR(ns);

			ncl = std::min(nsub, nclad);
			ncl_sqr = template_funcs::DSQR(ncl);

			aa = (nc_sqr / ns_sqr);

			bb = (nc_sqr / ncl_sqr);

			k = Two_PI / l;
			k_sqr = template_funcs::DSQR(k); // k_{0}^{2}

			k_sqr_nc_sqr = k_sqr * nc_sqr; // k_{0}^{2} n_{c}^{2}
			k_sqr_ns_sqr = k_sqr * ns_sqr; // k_{0}^{2} n_{s}^{2}
			k_sqr_ncl_sqr = k_sqr * ncl_sqr; // k_{0}^{2} n_{cl}^{2}

			double x = nc_sqr - ns_sqr;
			double y = ns_sqr - ncl_sqr;

			// WG asymmetry paramater
			g = y > 0.0 ? (y / x) : 0.0;

			na = sqrt(x); // numerical aperture

			V = (PI*d*na) / l; // V-parameter

			// predicted number of modes
			M = y > 0 ? static_cast<int>( std::max(1.0, ceil((2.0*V / PI) - atan(g) / PI)) ) : static_cast<int>(std::max(1.0, ceil((2.0*V / PI))));

			lower = k * ns; // lower bound of search space

			upper = k * nc; // upper bound of search space

			w = k * SPEED_OF_LIGHT;

			efieldint = 0.5*EPSILON*SPEED_OF_LIGHT;

			hfieldint = 0.5*MU*SPEED_OF_LIGHT;

			// Empty the std::vector each time a new instance of the class is called
			betaE.clear();
			betaH.clear();
		}
		else {
			std::string reason = "Error: void slab_tl_neff::set_params(double width,double lambda,double ncore,double nsub,double nclad)\n";
			if (!c1) reason += "WG width = " + template_funcs::toString(width, 3) + " is negative\n";
			if (!c2) reason += "Wavelength = " + template_funcs::toString(lambda, 3) + " is negative\n";
			if (!c3) reason += "Cladding Index = " + template_funcs::toString(nclad, 3) + " is less than one\n";
			if (!c4) reason += "Substrate Index = " + template_funcs::toString(nsub, 3) + " is less than one\n";
			if (!c5) reason += "Core Index = " + template_funcs::toString(ncore, 3) + " is less than cladding / substrate index\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl_neff::eigeneqn_3(double x, bool t, int mm)
{
	// Using this version of the dispersion equation means that you don't have to do a bracketing step
	// which you would have to do if you used eigeneqn

	try {
		
		if (k_sqr_nc_sqr > k_sqr_ns_sqr) {

			double h, p, q, tmp;

			double x_sqr = template_funcs::DSQR(x);

			tmp = k_sqr_nc_sqr - x_sqr;
			h = (tmp>0 ? sqrt(tmp) : 0.0);

			tmp = x_sqr - k_sqr_ns_sqr;
			p = (tmp>0 ? sqrt(tmp) : 0.0);

			tmp = x_sqr - k_sqr_ncl_sqr;
			q = (tmp>0 ? sqrt(tmp) : 0.0);

			if (t) {//TE modes
				return d * h - mm * PI - atan(q / h) - atan(p / h);
			}
			else {//TM modes
				return d * h - mm * PI - atan(bb*(q / h)) - atan(aa*(p / h));
			}
		}
		else {
			std::string reason = "Error: double slab_tl_neff::eigeneqn_3(double x,bool t,int mm)\n";
			reason += "Input parameters not correctly defined\n";
			reason += "k_{0}^{2} n_{c}^{2} = " + template_funcs::toString(k_sqr_nc_sqr) + ", k_{0}^{2} n_{s}^{2} = " + template_funcs::toString(k_sqr_ns_sqr) + "\n";
			return 0;
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl_neff::zbrent(double x1, double x2, double tol, bool t, int mm)
{
	//Find the root of the function between x1 and x2 using brent's method
	//The root is refined to +/- tol
	//Seemingly one of the best methods to use

	//This will be used to compute the roots of eigeneq_3
	//R. Sheehan 28 - 5 - 2010

	try {

		bool c1 = fabs(x1 - x2) > 1.0e-9 ? true : false; // cannot have x1 == x2
		bool c2 = tol > 1.0e-16 ? true : false;
		bool c3 = mm < M ? true : false;

		if (c1 && c2 && c3) {

			int iter;

			static const int ITMAX = 100;//Used in rtbis, rtflsp, rtsec, zriddr

			double a = std::min(x1, x2), b = std::max(x1, x2), c = std::max(x1, x2), d, e, min1, min2;
			double fc, p, q, r, s, tol1, xm;
			double fa = eigeneqn_3(a, t, mm), fb = eigeneqn_3(b, t, mm);

			if ((fa>0.0 && fb>0.0) || (fa<0.0 && fb<0.0)) {
				std::cerr << "Root must be bracketed in zbrent\n";
			}
			fc = fb;
			for (iter = 1; iter <= ITMAX; iter++) {
				if ((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0)) {
					c = a;
					fc = fa;
					e = d = b - a;
				}
				if (fabs(fc)<fabs(fb)) {
					a = b;
					b = c;
					c = a;
					fa = fb;
					fb = fc;
					fc = fa;
				}
				tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
				xm = 0.5*(c - b);
				if (fabs(xm) <= tol1 || fb == 0.0) return b;
				/*if(fabs(xm)<=tol1 || fb==0.0){
				std::cout<<"Brent's Method converged in "<<iter<<" iterations\n";
				return b;
				}*/
				if (fabs(e) >= tol1 && fabs(fa)>fabs(fb)) {
					s = fb / fa;
					if (a == c) {
						p = 2.0*xm*s;
						q = 1.0 - s;
					}
					else {
						q = fa / fc;
						r = fb / fc;
						p = s * (2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
						q = (q - 1.0)*(r - 1.0)*(s - 1.0);
					}
					if (p>0.0) q = -q;
					p = fabs(p);
					min1 = 3.0*xm*q - fabs(tol1*q);
					min2 = fabs(e*q);
					if (2.0*p<std::min(min1, min2)) {
						e = d;
						d = p / q;
					}
					else {
						d = xm;
						e = d;
					}
				}
				else {
					d = xm;
					e = d;
				}
				a = b;
				fa = fb;
				if (fabs(d)>tol1) {
					b += d;
				}
				else {
					b += template_funcs::SIGN(tol1, xm);
				}
				fb = eigeneqn_3(b, t, mm);
			}
			std::cerr << "Maximum number of iterations exceeded in zbrent\n";
			return 0.0;
		}
		else {
			std::string reason = "Error: double slab_tl_neff::zbrent(double x1,double x2,double tol,bool t,int mm)\n";
			if (!c1) reason += "Cannot have x1 = x2\nx1 = " + template_funcs::toString(x1) + ", x2 = " + template_funcs::toString(x1) + "\n";
			if (!c2) reason += "Desired tolerance is less than smallest allowable EPS\n";
			if (!c3) reason += "mm >= M\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_tl_neff::neff_search(bool mode)
{
	//This version solves the standard slab dispersion equation with a superior root finder
	//sic. Brent's Method NRinC ch. 9
	//R. Sheehan 28 - 5 - 2010

	try {

		if (lower < upper) {

			int m;
			double b;

			std::vector<double> vec;

			for (m = 0; m < M; m++) {
				b = zbrent(lower, upper, EPS, mode, m); 
				
				if (b>lower && b<upper) {
					vec.push_back(b);
				}

			}

			if (mode) {
				betaE = vec;
			}
			else {
				betaH = vec;
			}

			vec.clear();
		}
		else {
			std::string reason = "Error: void slab_tl_neff::neff_search(bool mode)\n";
			reason += "Search range is not correctly defined\n";
			reason += "lower = " + template_funcs::toString(lower, 4) + ", upper = " + template_funcs::toString(upper, 4) + "\n";
			throw std::range_error(reason);
		}
	}
	catch (std::range_error &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_tl_neff::report(bool mode)
{
	// print the computed results to screen
	std::cout << "Output for the slab waveguide calculator\n";
	std::cout << "\n";
	std::cout << "Waveguide width = " << d << " microns\n";
	std::cout << "Wavelength = " << l << " microns\n";
	std::cout << "Wavenumber = " << k << " inverse microns\n";
	std::cout << "\n";
	std::cout << "Core Index = " << nc << "\n";
	std::cout << "Substrate Index = " << ns << "\n";
	std::cout << "Cladding Index = " << ncl << "\n";
	std::cout << "\n";
	std::cout << "Numerical Aperture = " << na << "\n";
	std::cout << "V-Parameter = " << V << "\n";
	std::cout << "Number of Modes = " << M << "\n";
	std::cout << "\n";

	std::string pol = (mode ? "TE" : "TM");

	std::cout << pol << " Modes\n";
	std::cout << "There are " << nbeta(mode) << " calculated " + pol + " modes\n";

	for (int i = 0; i < nbeta(mode); i++) {
		std::cout << pol + "_{" << i + 1 << "} = " << std::setprecision(10) << _beta(i, mode) << " , n_eff_{" << i + 1 << "} = " << std::setprecision(10) << (_beta(i, mode) / k) << "\n";
	}
	std::cout << "\n";
}

int slab_tl_neff::get_nmodes(bool mode)
{
	// return the number of computed propagation constants for a given polarisation
	// it can be the case that the predicted value M is greater than the actual number of modes in a waveguide

	return nbeta(mode);
}

double slab_tl_neff::get_neff(int i, bool mode)
{
	// return the i^{th} effective index for a given polarisation
	// by convention mode == true => TE modes, mode == false => TM modes

	try {
		if (nbeta(mode) > 0) {
			if (i<0 || i > nbeta(mode)) {
				throw std::range_error("Error: double slab_tl_neff::get_neff(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (mode) {
					return betaE[i] / k;
				}
				else {
					return betaH[i] / k;
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl_neff::get_neff(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error & e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument & e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// Definitions for the three layer slab with optical mode profile calculation
slab_tl_mode::slab_tl_mode(void)
{
	// Default Constructor

}

slab_tl_mode::slab_tl_mode(double width, double lambda, double ncore, double nsub, double nclad)
{
	// Default Constructor
	set_params(width, lambda, ncore, nsub, nclad);
}

void slab_tl_mode::compute_neff(bool mode)
{
	// Compute the effective indices of the slab waveguide
	neff_search(mode); 

	std::string pol = (mode ? "TE" : "TM");

	std::cout << "There are " << nbeta(mode) << " calculated " + pol + " modes\n";
	for (int i = 0; i<static_cast<int>(nbeta(mode)); i++) {
		std::cout << "beta[" << i + 1 << "] = " << std::setprecision(6) << _beta(i, mode) << " , n_{eff} = " << _beta(i, mode) / k << "\n";
	}
	std::cout << "\n";
}

void slab_tl_mode::output_modes(bool mode, int N, double Lx, std::string &storage_directory)
{
	// This function will calculate the solutions corresponding to each mode and output them to a file for a single polarisation

	double xi;

	std::string filename;

	std::string pol = (mode ? "TE" : "TM");

	std::ofstream write;

	if (nbeta(mode) > 0) {

		std::vector< std::vector< double > > mat;

		mat.resize(nbeta(mode) + 2); // We want to output M solutions plus the corresponding positions

		double dx = Lx / (static_cast<double>(N - 1));

		xi = -0.5*Lx;

		mat[0].resize(N + 2);

		for (int j = 1; j <= N; j++) {
			mat[0][j] = xi;
			xi += dx;
		}

		for (int i = 1; i <= nbeta(mode); i++) {
			mat[i].resize(N + 2);
			for (int j = 1; j <= N; j++) {
				mat[i][j] = TE_TM(mat[0][j], i - 1, mode);
			}
		}

		// Output all the modes to a single file
		filename = storage_directory + pol + "_Mode_Profiles.txt";
		write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);

		if (write.is_open()) {

			for (int i = 1; i <= N; i++) {
				for (int j = 0; j<(nbeta(mode) + 1); j++)
					if (j == nbeta(mode))
						write << std::setprecision(20) << mat[j][i];
					else
						write << std::setprecision(20) << mat[j][i] << " , ";
				write << "\n";
			}

			write.close();

		}

		// output the sine-cosine form of the dispersion equation
		double db = ((upper - lower) / (99));

		filename = storage_directory + pol + "_Dispersion_Eqn.txt";
		write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);

		for (int i = 1; i<99; i++) {
			write << lower + i * db << " , " << std::setprecision(20) << eigeneqn(lower + i * db, mode) << "\n";
		}

		write.close();
	}
}

void slab_tl_mode::output_all_stats(std::string &storage_directory)
{
	//This function outputs all the stats associated with a particular calculation

	std::string file;

	file = storage_directory + "Slab_WG_Stats.txt";

	std::ofstream stats;
	stats.open(file.c_str(), std::ios_base::out | std::ios_base::trunc);

	stats << "Output File for the slab waveguide calculator\n";
	stats << "\n";
	stats << "Waveguide width = " << d << " microns\n";
	stats << "Wavelength = " << l << " microns\n";
	stats << "Wavenumber = " << k << " inverse microns\n";
	stats << "\n";
	stats << "Core Index = " << nc << "\n";
	stats << "Substrate Index = " << ns << "\n";
	stats << "Cladding Index = " << ncl << "\n";
	stats << "\n";
	stats << "Numerical Aperture = " << na << "\n";
	stats << "V-Parameter = " << V << "\n";
	stats << "Number of Modes = " << M << "\n";
	stats << "\n";

	if (nbeta(TE) > 0) {

		output_stats(TE, stats);

	}

	if (nbeta(TM) > 0) {

		output_stats(TM, stats);

	}

	stats.close();
}

double slab_tl_mode::g1(int i, bool t)
{
	try {
		
		if (nbeta(t)) {

			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl_mode::g1(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					double beta_sqr = template_funcs::DSQR(betaE[i]);
					return (beta_sqr / k_sqr_nc_sqr) + (beta_sqr / k_sqr_ns_sqr) - 1.0;
				}
				else {
					double beta_sqr = template_funcs::DSQR(betaH[i]);
					return (beta_sqr / k_sqr_nc_sqr) + (beta_sqr / k_sqr_ns_sqr) - 1.0;
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl_mode::g1(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl_mode::g2(int i, bool t)
{
	try {

		if (nbeta(t) > 0) {

			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl_mode::g2(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					double beta_sqr = template_funcs::DSQR(betaE[i]);
					return (beta_sqr / k_sqr_nc_sqr) + (beta_sqr / k_sqr_ncl_sqr) - 1.0;
				}
				else {
					double beta_sqr = template_funcs::DSQR(betaH[i]);
					return (beta_sqr / k_sqr_nc_sqr) + (beta_sqr / k_sqr_ncl_sqr) - 1.0;
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl_mode::g2(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl_mode::deff(int i, bool t)
{
	// Effective width of waveguide mode

	try {

		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl_mode::deff(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return d + (1.0 / p(i, t)) + (1.0 / q(i, t));
				}
				else {
					return d + (1.0 / (g1(i, t) * p(i, t))) + (1.0 / (g2(i, t) * q(i, t)));
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl_mode::deff(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl_mode::phase(int i, bool t)
{
	// Phase of mode in slab waveguide

	try {
		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl_neff::phase(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					double hh = h(i, t);
					double pp = p(i, t);
					double qq = q(i, t);
					return i * (PI / 2) + 0.5 * atan((pp / hh)) - 0.5 * atan((qq / hh));
				}
				else {
					double hh = h(i, t);
					double pp = p(i, t);
					double qq = q(i, t);
					return  i * (PI / 2) + 0.5*atan(aa*(pp / hh)) - 0.5*atan(bb * (qq / hh));
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl_mode::phase(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl_mode::norm_const(int i, bool t)
{
	// Normalisation constant for mode in a slab waveguide

	try {

		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl_mode::norm_const(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt((w * MU) / (betaE[i] * deff(i, t)));
				}
				else {
					return sqrt((w * EPSILON * nc_sqr) / (betaH[i] * deff(i, t)));
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl_mode::norm_const(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl_mode::conf_fact(int i, bool t)
{
	// Confinement factor of mode in slab waveguide

	try {

		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl_neff::conf_fact(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					double hh_sqr = template_funcs::DSQR(h(i, t));
					double pp = p(i, t);
					double pp_sqr = template_funcs::DSQR(pp);
					double qq = q(i, t);
					double qq_sqr = template_funcs::DSQR(qq);

					return (d + (1.0 / pp)*(1.0 / (1.0 + hh_sqr / pp_sqr)) + (1.0 / qq)*(1.0 / (1.0 + hh_sqr / qq_sqr))) / (deff(i, t));
				}
				else {
					double hh_sqr = template_funcs::DSQR(h(i, t));
					double pp = p(i, t);
					double pp_sqr = template_funcs::DSQR(pp);
					double qq = q(i, t);
					double qq_sqr = template_funcs::DSQR(qq);
					return (d + (1.0 / (g1(i, t) * pp))*(1.0 / (1.0 + hh_sqr / pp_sqr)) + (1.0 / (g2(i, t) * qq))*(1.0 / (1.0 + hh_sqr / qq_sqr))) / (deff(i, t));
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl_mode::conf_fact(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl_mode::eigeneqn(double x, bool t)
{
	// Dispersion equation for slab waveguide written as a sum of sine and cosine functions
	// The version in eigeneqn_3, written as a sum over atan(*) is preferred for fast root finding

	try {

		if (x >= lower && x <= upper) {

			double x_sqr = template_funcs::DSQR(x);
			double h = sqrt(k_sqr_nc_sqr - x_sqr);
			double p = sqrt(x_sqr - k_sqr_ns_sqr);
			double q = sqrt(x_sqr - k_sqr_ncl_sqr);

			if (t) {//TE modes
				return (template_funcs::DSQR(h) - (p*q)) * sin((d * h)) - h * (p + q) * cos((d * h));
			}
			else {//TM modes
				return (template_funcs::DSQR(h) - aa * bb * (p * q)) * sin((d * h)) - (h * aa * p + h * bb * q) * cos((d * h));
			}
		}
		else {
			return 0;
			throw std::range_error("Error: double slab_tl_mode::eigeneqn(double x,bool t)\n Attempting to compute dispersion equation outside valid range\n");
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
}

double slab_tl_mode::TE_TM(double x, int i, bool mode)
{
	// Function which defines the shape of the modes
	// This method uses the stored computed propagation constants
	// Assumes that the slab wg has core region defined on -t < x < t, t = d/2

	try {

		if (nbeta(mode) > 0) {

			if (i<0 || i > nbeta(mode) ) {
				throw std::range_error("Error: double slab_tl_mode::TE_TM(double x,int i,bool mode)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {

				if (x > tt) {
					return norm_const(i, mode)*cos(h(i, mode)*tt - phase(i, mode))*exp(-q(i, mode)*(x - tt)); // Cladding
				}
				else if (x < -tt) {
					return norm_const(i, mode)*cos(h(i, mode)*tt + phase(i, mode))*exp(p(i, mode)*(x + tt)); // Substrate
				}
				else {
					return norm_const(i, mode)*cos(h(i, mode)*x - phase(i, mode)); // Core
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl_mode::TE_TM(int i, bool mode)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_tl_mode::output_stats(bool mode, std::ofstream &file_obj)
{
	// output the results for a particular polarisation
	// R. Sheehan 18 - 7 - 2014

	std::string pol = (mode ? "TE" : "TM");

	if (file_obj.is_open()) {

		file_obj << pol << " Modes\n";
		file_obj << "There are " << nbeta(mode) << " calculated " + pol + " modes\n";

		for (int i = 0; i<nbeta(mode); i++) {
			file_obj << pol + "_{" << i + 1 << "} = " << std::setprecision(10) << _beta(i, mode) << " , n_eff_{" << i + 1 << "} = " << std::setprecision(10) << (_beta(i, mode) / k) << "\n";
		}
		file_obj << "\n";
		for (int i = 0; i<nbeta(mode); i++) {
			file_obj << "Confinement Factor " << i + 1 << " = " << std::setprecision(10) << conf_fact(i, mode) << "\n";
		}
		file_obj << "\n";
		for (int i = 0; i<nbeta(mode); i++) {
			file_obj << "Normalisation Constant " << i + 1 << " = " << std::setprecision(10) << norm_const(i, mode) << "\n";
		}
		file_obj << "\n";
		for (int i = 0; i<nbeta(mode); i++) {
			file_obj << "Phase " << i + 1 << " = " << std::setprecision(10) << phase(i, mode) << "\n";
		}
		file_obj << "\n";
		for (int i = 0; i<nbeta(mode); i++) {
			file_obj << "Effective Width " << i + 1 << " = " << std::setprecision(10) << deff(i, mode) << "\n";
		}
		file_obj << "\n";
		for (int i = 0; i<nbeta(mode); i++) {
			file_obj << "h_" << i + 1 << " = " << std::setprecision(10) << h(i, mode) << "\n";
		}
		file_obj << "\n";
		for (int i = 0; i<nbeta(mode); i++) {
			file_obj << "p_" << i + 1 << " = " << std::setprecision(10) << p(i, mode) << "\n";
		}
		file_obj << "\n";
		for (int i = 0; i<nbeta(mode); i++) {
			file_obj << "q_" << i + 1 << " = " << std::setprecision(10) << q(i, mode) << "\n";
		}
		file_obj << "\n";
	}
}

//double slab_tl_neff::coupling_coefficient(double &separ, bool mode)
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
//		if (nbeta(mode) > 0) {
//
//			double tt = d / 2;
//			double hh = h(0, mode);
//			double pp = p(0, mode);
//			double vs = template_funcs::DSQR(k * na);
//			double num = template_funcs::DSQR(hh * pp);
//			double denom = _beta(0, mode)*(1.0 + (tt * pp))*vs;
//			double arg = -1.0 * pp * (separ - d);
//
//			if (denom > 0) {
//				return ((num / denom)* exp(arg));
//			}
//			else {
//				std::string reason = "Error: double slab_tl_neff::coupling_coefficient(double &separ, bool mode)\n"; 
//				reason += "Division by zero occurred\n"; 
//				throw std::runtime_error(reason);
//				return 0.0;
//			}
//		}
//		else {
//			std::string reason = "Error: double slab_tl_neff::coupling_coefficient(double &separ, bool mode)\n"; 
//			reason += "No modes computed for this waveguide\n";
//			throw std::runtime_error(reason);
//			return 0.0;
//		}
//	}
//	catch (std::runtime_error &e) {
//		useful_funcs::exit_failure_output(e.what());
//		exit(EXIT_FAILURE);
//	}
//}

// Definitions for the four layer slab derived class

slab_fl_neff_A::slab_fl_neff_A()
{
	// Default Constructor
}

slab_fl_neff_A::slab_fl_neff_A(double width, double rib_width, double lambda, double ncore, double nsub, double nclad, double nrib)
{
	// Constructor

	// Case A => Field Oscillating in Core and Ridge
	// For there to be a solution one has to have ns <= ncl < nr < nc

	set_params(width, rib_width, lambda, ncore, nsub, nclad, nrib);
}

void slab_fl_neff_A::set_params(double width, double rib_width, double lambda, double ncore, double nsub, double nclad, double nrib)
{
	// set parameters for four layer slab
	
	// Case A => Field Oscillating in Core and Ridge
	// For there to be a solution one has to have ns <= ncl < nr < nc

	// Case B: Field Oscillating in Core Only
	// For there to be a solution one has to have ncl < nm < nc, where nm = Max(nr,ns)

	try {
		bool c1 = width > 0.0 ? true : false;
		bool c1a = rib_width > 0.0 ? true : false; 
		bool c2 = lambda > 0.0 ? true : false;
		bool c3 = nclad >= 1.0 ? true : false;
		bool c4 = nsub >= 1.0 ? true : false;
		bool c6 = nrib >= 1.0 ? true : false;
		bool c5 = ncore > std::max(nsub, nclad) ? true : false;
		bool c7 = ncore > std::max(nsub, nrib) ? true : false;

		if (c1 && c1a && c2 && c3 && c4 && c5 && c6 && c7) {

			d = width;

			tt = 0.5 * d;

			dr = rib_width; 

			l = lambda;

			nc = ncore;
			nc_sqr = template_funcs::DSQR(nc);

			nr = nrib;
			nr_sqr = template_funcs::DSQR(nr);

			ns = nsub;
			ns_sqr = template_funcs::DSQR(ns);

			ncl = nclad;
			ncl_sqr = template_funcs::DSQR(ncl);

			nm = std::max(ns, ncl);
			nm_sqr = template_funcs::DSQR(nm);

			k = Two_PI / l;
			k_sqr = template_funcs::DSQR(k); // k_{0}^{2}

			k_sqr_nc_sqr = k_sqr * nc_sqr; // k_{0}^{2} n_{c}^{2}
			k_sqr_ns_sqr = k_sqr * ns_sqr; // k_{0}^{2} n_{s}^{2}
			k_sqr_ncl_sqr = k_sqr * ncl_sqr; // k_{0}^{2} n_{cl}^{2}
			k_sqr_nr_sqr = k_sqr * nr_sqr; // k_{0}^{2} n_{r}^{2}
			k_sqr_nm_sqr = k_sqr * nm_sqr; // k_{0}^{2} n_{m}^{2}

			etacs = nc_sqr / ns_sqr;
			etacr = nc_sqr / nr_sqr;
			etarcl = nr_sqr / ncl_sqr;

			// Only difference between A, B cases is search space for neff and eigenequation
			// Case A: lower = k max{ncl, ns}, upper = k nr, NA^{2} = nr^{2} - max{ncl, ns}^{2}
			// Case B: lower = k nm, upper =k nc, NA^{2} = nc^{2} - nm^{2}

			double x = nr_sqr - nm_sqr;

			na = sqrt(x); // numerical aperture

			V = (PI * (dr) * na) / l; // V-parameter

			// predicted number of modes
			M = static_cast<int>( std::max( 1.0, ceil( (2.0*V / PI) ) ) );

			lower = k * nm; // lower bound of search space for case A

			upper = k * nr; // upper bound of search space for case A

			w = k * SPEED_OF_LIGHT;

			efieldint = 0.5*EPSILON*SPEED_OF_LIGHT;

			hfieldint = 0.5*MU*SPEED_OF_LIGHT;

			// Empty the std::vector each time a new instance of the class is called
			betaE.clear();
			betaH.clear();
		}
		else {
			std::string reason = "Error: void slab_fl_neff_A::set_params(double width,double lambda,double ncore,double nsub,double nclad)\n";
			if (!c1) reason += "WG width = " + template_funcs::toString(width, 3) + " is negative\n";
			if (!c2) reason += "Wavelength = " + template_funcs::toString(lambda, 3) + " is negative\n";
			if (!c3) reason += "Cladding Index = " + template_funcs::toString(nclad, 3) + " is less than one\n";
			if (!c4) reason += "Substrate Index = " + template_funcs::toString(nsub, 3) + " is less than one\n";
			if (!c6) reason += "Rib Index = " + template_funcs::toString(nrib, 3) + " is less than one\n";
			if (!c5) reason += "Core Index = " + template_funcs::toString(ncore, 3) + " is less than cladding / substrate index\n";
			if (!c7) reason += "Core Index = " + template_funcs::toString(ncore, 3) + " is less than rib / substrate index\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_fl_neff_A::eigeneqn_a(double x, int mm, bool t)
{
	// Dispersion equation corresponding to case a from Adams
	
	// Case A => Field Oscillating in Core and Ridge
	// For there to be a solution one has to have ns <= ncl < nr < nc

	// https://xkcd.com/2034/

	try{
		
		if (k_sqr_nr_sqr > k_sqr_nm_sqr) {

			double h, p, q, r, tmp;

			double x_sqr = template_funcs::DSQR(x);

			tmp = k_sqr_nc_sqr - x_sqr;
			h = (tmp > 0 ? sqrt(tmp) : 0.0);

			tmp = k_sqr_nr_sqr - x_sqr;
			r = (tmp>0 ? sqrt(tmp) : 0.0);

			tmp = x_sqr - k_sqr_ns_sqr;
			p = (tmp>0 ? sqrt(tmp) : 0.0);

			tmp = x_sqr - k_sqr_ncl_sqr;
			q = (tmp>0 ? sqrt(tmp) : 0.0);

			if (t) {//TE modes
				return (d * h) - atan( (p / h) ) - atan( (r / h) * tan( atan(q / r) - (dr * r) ) ) - mm * PI;
			}
			else {//TM modes
				return (d * h) - atan( etacs*(p / h) ) - atan( etacr * (r / h) * tan( atan( etarcl*(q / r) ) - (dr * r) ) ) - mm * PI;
			}
		}
		else {
			std::string reason = "Error: double slab_fl_neff_A::eigeneqn_a(double x, int mm, bool t)\n";
			reason += "Input parameters not correctly defined\n";
			reason += "k_{0}^{2} n_{c}^{2} = " + template_funcs::toString(k_sqr_nc_sqr) + ", k_{0}^{2} n_{r}^{2} = " + template_funcs::toString(k_sqr_nr_sqr) + "\n";
			return 0;
		}	
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_fl_neff_A::eigeneqn_aa(double x, bool t)
{
	// Dispersion equation corresponding to case a from Adams

	// Case A => Field Oscillating in Core and Ridge
	// For there to be a solution one has to have ns <= ncl < nr < nc

	// equation re-written in an attempt to ameliorate effects of discontinuities arising from taking tangent of various items
	// R. Sheehan 18 - 12 - 2020

	try {

		if (k_sqr_nr_sqr > k_sqr_nm_sqr) {

			double h, p, q, r, tmp, Wh, psi, t1, t2, t3, t4;

			double x_sqr = template_funcs::DSQR(x);

			tmp = k_sqr_nc_sqr - x_sqr;
			h = (tmp > 0 ? sqrt(tmp) : 0.0);

			tmp = k_sqr_nr_sqr - x_sqr;
			r = (tmp > 0 ? sqrt(tmp) : 0.0);

			tmp = x_sqr - k_sqr_ns_sqr;
			p = (tmp > 0 ? sqrt(tmp) : 0.0);

			tmp = x_sqr - k_sqr_ncl_sqr;
			q = (tmp > 0 ? sqrt(tmp) : 0.0);

			Wh = d * h;
			t1 = sin(Wh);
			t2 = cos(Wh);

			if (t) {//TE modes
				t3 = (p / h);
				t4 = (r / h);
				psi = atan(q / r) - (dr * r);
			}
			else {//TM modes
				t3 = etacs * (p / h);
				t4 = etacr * (r / h);
				psi = atan( etarcl * (q / r) ) - (dr * r);
			}

			return ( ( t1 - (t3 * t2) ) - ( (t3 * t1) + t2) * t4 * tan(psi) );
		}
		else {
			std::string reason = "Error: double slab_fl_neff_A::eigeneqn_a(double x, int mm, bool t)\n";
			reason += "Input parameters not correctly defined\n";
			reason += "k_{0}^{2} n_{c}^{2} = " + template_funcs::toString(k_sqr_nc_sqr) + ", k_{0}^{2} n_{r}^{2} = " + template_funcs::toString(k_sqr_nr_sqr) + "\n";
			return 0;
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_fl_neff_A::zbrent(double x1, double x2, double tol, bool t, int mm)
{
	//Find the root of the function between x1 and x2 using brent's method
	//The root is refined to +/- tol
	//Seemingly one of the best methods to use

	//This will be used to compute the roots of eigeneq_b
	//R. Sheehan 28 - 5 - 2010

	try {

		bool c1 = fabs(x1 - x2) > 1.0e-9 ? true : false; // cannot have x1 == x2
		bool c2 = tol > 1.0e-16 ? true : false;
		bool c3 = mm < M ? true : false;

		if (c1 && c2 && c3) {

			int iter;

			static const int ITMAX = 100;//Used in rtbis, rtflsp, rtsec, zriddr

			double a = std::min(x1, x2), b = std::max(x1, x2), c = std::max(x1, x2), d, e, min1, min2;
			double fc, p, q, r, s, tol1, xm;
			//double fa = eigeneqn_a(a, mm, t), fb = eigeneqn_a(b, mm, t);
			double fa = eigeneqn_aa(a, t), fb = eigeneqn_aa(b, t);

			if ((fa>0.0 && fb>0.0) || (fa<0.0 && fb<0.0)) {
				std::cerr << "Root must be bracketed in zbrent\n";
			}
			fc = fb;
			for (iter = 1; iter <= ITMAX; iter++) {
				if ((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0)) {
					c = a;
					fc = fa;
					e = d = b - a;
				}
				if (fabs(fc)<fabs(fb)) {
					a = b;
					b = c;
					c = a;
					fa = fb;
					fb = fc;
					fc = fa;
				}
				tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
				xm = 0.5*(c - b);
				if (fabs(xm) <= tol1 || fb == 0.0) return b;
				if(fabs(xm)<=tol1 || fb==0.0){
					std::cout<<"Brent's Method converged in "<<iter<<" iterations\n";
					return b;
				}
				if (fabs(e) >= tol1 && fabs(fa)>fabs(fb)) {
					s = fb / fa;
					if (a == c) {
						p = 2.0*xm*s;
						q = 1.0 - s;
					}
					else {
						q = fa / fc;
						r = fb / fc;
						p = s * (2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
						q = (q - 1.0)*(r - 1.0)*(s - 1.0);
					}
					if (p>0.0) q = -q;
					p = fabs(p);
					min1 = 3.0*xm*q - fabs(tol1*q);
					min2 = fabs(e*q);
					if (2.0*p<std::min(min1, min2)) {
						e = d;
						d = p / q;
					}
					else {
						d = xm;
						e = d;
					}
				}
				else {
					d = xm;
					e = d;
				}
				a = b;
				fa = fb;
				if (fabs(d)>tol1) {
					b += d;
				}
				else {
					b += template_funcs::SIGN(tol1, xm);
				}
				//fb = eigeneqn_a(b, mm, t);
				fb = eigeneqn_aa(b, t);
			}
			std::cerr << "Maximum number of iterations exceeded in zbrent\n";
			return 0.0;
		}
		else {
			std::string reason = "Error: double slab_fl_neff_A::zbrent(double x1,double x2,double tol,bool t,int mm)\n";
			if (!c1) reason += "Cannot have x1 = x2\nx1 = " + template_funcs::toString(x1) + ", x2 = " + template_funcs::toString(x1) + "\n";
			if (!c2) reason += "Desired tolerance is less than smallest allowable EPS\n";
			if (!c3) reason += "mm >= M\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

//void slab_fl_neff_A::neff_search(bool mode)
//{
//	// Compute the waveguide mode effective indices based on given polarisation and wg_type
//	
//	// Case A => Field Oscillating in Core and Ridge
//	// For there to be a solution one has to have ns <= ncl < neff < nr
//
//	// Case A: lower = k ncl, upper = k nr, NA^{2} = nr^{2} - ncl^{2}
//
//	try {
//		if (lower < upper) {
//		
//			int m;
//			double b;
//
//			std::vector<double> vec;
//
//			for (m = 0; m < M; m++) {
//				b = zbrent(lower, upper, EPS, mode, m);
//
//				if (b>lower && b<upper) {
//					vec.push_back(b);
//				}
//
//			}
//
//			if (mode) {
//				betaE = vec;
//			}
//			else {
//				betaH = vec;
//			}
//
//			vec.clear();
//		}
//		else {
//			std::string reason = "Error: void slab_fl_neff_A::neff_search(bool mode, bool wg_type)\n";
//			reason += "Search range is not correctly defined\n";
//			reason += "lower = " + template_funcs::toString(lower, 4) + ", upper = " + template_funcs::toString(upper, 4) + "\n";
//			throw std::range_error(reason);
//		}	
//	}
//	catch (std::range_error &e) {
//		useful_funcs::exit_failure_output(e.what());
//		exit(EXIT_FAILURE);
//	}
//}

void slab_fl_neff_A::neff_search(bool mode) 
{
	// Alternative neff_search based on fuinding roots of eigeneqn_aa
	// R. Sheehan 4 - 1 - 2021

	// Case A => Field Oscillating in Core and Ridge
	// For there to be a solution one has to have nm < neff < nr

	// Case A: lower = k nm, upper = k nr, NA^{2} = nr^{2} - nm^{2}

	// Roots must be bracket before proceeding

	try {
		bool c1 = lower < upper ? true : false; 
		bool c2 = brackets.size() > 0 ? true : false; 
		bool c10 = c1 && c2; 
		if (c10) {
			int m = 0;
			double b;

			std::vector<double> vec;

			for (int i = static_cast<int>( brackets.size() ) - 1; i >= 0 ; i--) {

				b = zbrent(brackets[i].get_x_lower(), brackets[i].get_x_upper(), EPS, mode, m);

				if (b > lower && b < upper) {
					vec.push_back(b);
				}
			}

			if (mode) {
				betaE = vec;
			}
			else {
				betaH = vec;
			}

			vec.clear();
		}
		else {
			std::string reason = "Error: void slab_fl_neff_A::neff_search_aa(bool mode, bool wg_type)\n";
			if(!c1) reason += "Search range is not correctly defined\n";
			if(!c1) reason += "lower = " + template_funcs::toString(lower, 4) + ", upper = " + template_funcs::toString(upper, 4) + "\n";
			if (!c2) reason += "Bracketing intervals for roots have not been found\n"; 
			throw std::range_error(reason);
		}
	}
	catch (std::range_error& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_fl_neff_A::bracket_roots(bool mode, bool loud)
{
	// step through the search space and see what intervals contain roots of the dispersion equation
	// Bracketing algorithm taken from earlier work https://github.com/rnsheehan/find_roots_1D
	// R. Sheehan 4 - 1 - 2021

	try {
		if (lower < upper) {
		
			int N = 201; // divide the search space into a number of sub-intervals

			// compute the step-size in the search space
			double dx = (upper - lower) / (static_cast<double>(N - 1));

			double xl = lower + dx;
			double xu = xl + dx; 
			double fl, fu; 

			// remove any stored intervals
			brackets.clear(); 

			for (int i = 1; i <= N - 2; i++) {

				// Evaluate the dispersion equation at the interval endpoints
				if (i == 1) {
					fl = template_funcs::Signum( eigeneqn_aa(xl, mode) );
				}
				else {
					fl = fu; 
				}

				fu = template_funcs::Signum( eigeneqn_aa(xu, mode) );

				// Perform the bisection test
				if (fl * fu < 0.0) {
					// the sub-interval contains a root so it is stored
					brackets.push_back(interval(xl, xu));
				}

				// update the endpoints of the sub-interval
				xl = xu;
				xu += dx; 
			}

			if (loud && brackets.size() > 0) {
				std::cout << "The search space contains " << brackets.size() << " roots\n"; 
				for (size_t j = 0; j < brackets.size(); j++) {
					std::cout << brackets[j].get_x_lower() << " , " << brackets[j].get_x_upper() << "\n"; 
				}
			}
		}
		else {
			std::string reason = "Error: void slab_fl_neff_A::bracket_roots(bool mode)\n";
			reason += "Search range is not correctly defined\n";
			reason += "lower = " + template_funcs::toString(lower, 4) + ", upper = " + template_funcs::toString(upper, 4) + "\n";
			throw std::range_error(reason);
		}
	}
	catch (std::range_error& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_fl_neff_A::output_disp_eqn_a(bool mode, std::string& storage_directory)
{
	// Output the case A four-layer slab dispersion equation
	// R. Sheehan 18 - 12 - 2020

	try {
		if (lower < upper) {
			int N = 101; 

			// create array to hold dispersion equations for all the modes
			std::vector< std::vector< double > > mat;

			mat.resize(M + 2); // We want to output M solutions plus the corresponding positions

			// store the positions at which the eigeneqn_a will be evaluated
			double dx = (upper - lower) / (static_cast<double>(N - 1));

			double xi = lower;

			mat[0].resize(N + 2);

			for (int j = 1; j <= N; j++) {
				mat[0][j] = xi;
				xi += dx;
			}

			// compute the values of the dispersion equation for each of the expected modes
			for (int i = 1; i <= M; i++) {
				mat[i].resize(N + 2);
				for (int j = 1; j <= N; j++) {
					mat[i][j] = eigeneqn_a(mat[0][j], i - 1, mode);
				}
			}

			// Output the computed dispersion equations
			std::string pol = (mode ? "TE" : "TM");
			std::string filename = storage_directory + pol + "_Disp_Eqns_FL_A.txt";
			std::ofstream write;
			
			write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);

			if (write.is_open()) {

				for (int i = 1; i <= N; i++) {
					for (int j = 0; j < (M + 1); j++)
						if (j == M)
							write << std::setprecision(20) << mat[j][i];
						else
							write << std::setprecision(20) << mat[j][i] << " , ";
					write << "\n";
				}

				write.close();
			}

			mat.clear(); 
		}
		else {
			std::string reason = "Error: void slab_fl_neff_A::output_disp_eqn(bool mode)\n";
			reason += "Search range is not correctly defined\n";
			reason += "lower = " + template_funcs::toString(lower, 4) + ", upper = " + template_funcs::toString(upper, 4) + "\n";
			throw std::range_error(reason);
		}
	}
	catch (std::range_error& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_fl_neff_A::output_disp_eqn_aa(bool mode, std::string& storage_directory)
{
	// Output the case A four-layer slab dispersion equation
	// R. Sheehan 18 - 12 - 2020

	try {
		if (lower < upper) {
			int N = 201;

			// store the positions at which the eigeneqn_a will be evaluated
			double dx = (upper - lower) / (static_cast<double>(N - 1));

			double xi = lower + dx;

			// Output the computed dispersion equations
			std::string pol = (mode ? "TE" : "TM");
			std::string filename = storage_directory + pol + "_Disp_Eqns_FL_AA.txt";
			std::ofstream write;

			write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);

			if (write.is_open()) {

				for (int i = 1; i <= N - 2; i++) {
					write << std::setprecision(10) << xi << " , " << eigeneqn_aa(xi, mode) << "\n"; 
					xi += dx; 
				}

				write.close();
			}
		}
		else {
			std::string reason = "Error: void slab_fl_neff_A::output_disp_eqn(bool mode)\n";
			reason += "Search range is not correctly defined\n";
			reason += "lower = " + template_funcs::toString(lower, 4) + ", upper = " + template_funcs::toString(upper, 4) + "\n";
			throw std::range_error(reason);
		}
	}
	catch (std::range_error& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// slab_fl_mode_A is used to compute the mode profile in the four layer slab case A

slab_fl_mode_A::slab_fl_mode_A()
{
	// Default
}

slab_fl_mode_A::slab_fl_mode_A(double width, double rib_width, double lambda, double ncore, double nsub, double nclad, double nrib)
{
	// Assign value to all the parameters
	set_params(width, rib_width, lambda, ncore, nsub, nclad, nrib);
}

void slab_fl_mode_A::compute_neff(bool mode)
{
	// Compute the effective indices of the slab waveguide
	neff_search(mode);

	std::string pol = (mode ? "TE" : "TM");

	std::cout << "Search Space: \n";
	std::cout << "lower: " << lower << " < beta < upper: " << upper << "\n"; 
	std::cout << "lower: " << nm << " < neff < upper: " << nr << "\n"; 
	std::cout << "There are " << M << " predicted modes\n"; 
	std::cout << "There are " << nbeta(mode) << " calculated " + pol + " modes\n";
	for (int i = 0; i<static_cast<int>(nbeta(mode)); i++) {
		std::cout << "beta[" << i + 1 << "] = " << std::setprecision(6) << _beta(i, mode) << " , n_{eff} = " << _beta(i, mode) / k << "\n";
	}
	std::cout << "\n";
}

double slab_fl_mode_A::phase(int i, bool t)
{
	// Phase of mode in type A four layer slab waveguide

	try {
		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_fl_mode_A::phase(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					double qq = q(i, t);
					double rr = rA(i, t);
					return ( atan( qq / rr ) - (dr * rr) );
				}
				else {
					double qq = q(i, t);
					double rr = rA(i, t);
					return ( atan( etarcl * ( qq / rr) ) - ( dr * rr ) );
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_fl_mode_B::phase(int i, bool t)s\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_fl_mode_A::TE_TM(double x, int i, bool mode)
{
	// shape of waveguide mode
	// Function which defines the shape of the modes for Type A four layer slab
	// This method uses the stored computed propagation constants
	// Assumes that the slab wg has core region defined on -W < x < 0 and 0 < x < D, 

	try {
		if (nbeta(mode) > 0) {
			if (i<0 || i > nbeta(mode)) {
				throw std::range_error("Error: double slab_fl_mode_A::TE_TM(double x,int i,bool mode)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				double hh = h(i, mode);
				double rr = rA(i, mode);
				double wh = d * hh;
				double Dr = dr * rr; 
				double qq = q(i, mode);
				double pp = p(i, mode);
				double ph = phase(i, mode);
				double rrhh = mode ? (rr / hh) * tan(ph) : (rr / hh) * tan(ph) * etacr;
				if (x < -d) {
					// x < -W
					return ((cos(wh) + (rrhh * sin(wh)))*exp(pp * (x + d)));
				}
				else if (x >= -d && x <= 0) {
					// -W < x < 0
					return (cos(hh * x) - rrhh * sin(hh * x));
				}
				else if (x > 0 && x <= dr) {
					// 0 < x < D
					return (cos((rr * x) + ph) / cos(ph));
				}
				else {
					// x > D
					return ((cos( Dr + ph) / cos(ph) )*exp(qq * (dr - x)));
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_fl_mode_A::TE_TM(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_fl_mode_A::output_modes(bool mode, int N, double Lx, std::string& storage_directory)
{
	// This function will calculate the solutions corresponding to each mode and output them to a file for a single polarisation

	double xi;

	std::string filename;

	std::string pol = (mode ? "TE" : "TM");

	std::ofstream write;

	if (nbeta(mode) > 0) {

		std::vector< std::vector< double > > mat;

		mat.resize(nbeta(mode) + 2); // We want to output M solutions plus the corresponding positions

		double dx = Lx / (static_cast<double>(N - 1));

		xi = -0.5 * Lx;

		mat[0].resize(N + 2);

		// store the position values in the array
		for (int j = 1; j <= N; j++) {
			mat[0][j] = xi;
			xi += dx;
		}

		// store the computed field values in the array
		for (int i = 1; i <= nbeta(mode); i++) {
			mat[i].resize(N + 2);
			for (int j = 1; j <= N; j++) {
				mat[i][j] = TE_TM(mat[0][j], i - 1, mode);
			}
		}

		// normalise the field values to have max-val = 1
		for (int i = 1; i <= nbeta(mode); i++) {
			// determine the norm of the mode stored in row i
			double norm = vecut::inf_norm(mat[i]); 
			for (int j = 1; j <= N; j++) {
				mat[i][j] /= norm;
			}
		}

		// Output all the modes to a single file
		filename = storage_directory + pol + "_Mode_Profiles_FL_A.txt";
		write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);

		if (write.is_open()) {

			for (int i = 1; i <= N; i++) {
				for (int j = 0; j < (nbeta(mode) + 1); j++)
					if (j == nbeta(mode))
						write << std::setprecision(20) << mat[j][i];
					else
						write << std::setprecision(20) << mat[j][i] << " , ";
				write << "\n";
			}

			write.close();

		}

		// output the dispersion equation for the case A four layer slab
		output_disp_eqn_aa(mode, storage_directory); 
	}
}

// slab_fl_mode_B is used to compute the effective indices in the four layer slab case B

slab_fl_neff_B::slab_fl_neff_B()
{
	// Default Constructor
}

slab_fl_neff_B::slab_fl_neff_B(double width, double rib_width, double lambda, double ncore, double nsub, double nclad, double nrib)
{
	// Constructor

	// Case B: Field Oscillating in Core Only
	// For there to be a solution one has to have ncl < nm < nc, where nm = Max(nr,ns)

	set_params(width, rib_width, lambda, ncore, nsub, nclad, nrib);
}

void slab_fl_neff_B::set_params(double width, double rib_width, double lambda, double ncore, double nsub, double nclad, double nrib)
{
	// set parameters for four layer slab

	// Case A => Field Oscillating in Core and Ridge
	// For there to be a solution one has to have ns <= ncl < nr < nc

	// Case B: Field Oscillating in Core Only
	// For there to be a solution one has to have ncl < nm < nc, where nm = Max(nr,ns)

	try {
		bool c1 = width > 0.0 ? true : false;
		bool c1a = rib_width > 0.0 ? true : false;
		bool c2 = lambda > 0.0 ? true : false;
		bool c3 = nclad >= 1.0 ? true : false;
		bool c4 = nsub >= 1.0 ? true : false;
		bool c6 = nrib >= 1.0 ? true : false;
		bool c5 = ncore > std::max(nsub, nclad) ? true : false;
		bool c7 = ncore > std::max(nsub, nrib) ? true : false;

		if (c1 && c1a && c2 && c3 && c4 && c5 && c6 && c7) {

			d = width;

			tt = 0.5 * d; 

			dr = rib_width;

			l = lambda;

			nc = ncore;
			nc_sqr = template_funcs::DSQR(nc);

			nr = nrib;
			nr_sqr = template_funcs::DSQR(nr);

			ns = std::max(nsub, nclad);
			ns_sqr = template_funcs::DSQR(ns);

			ncl = std::min(nsub, nclad);
			ncl_sqr = template_funcs::DSQR(ncl);

			nm = std::max(ns, nr);
			nm_sqr = template_funcs::DSQR(nm);

			k = Two_PI / l;
			k_sqr = template_funcs::DSQR(k); // k_{0}^{2}

			k_sqr_nc_sqr = k_sqr * nc_sqr; // k_{0}^{2} n_{c}^{2}
			k_sqr_ns_sqr = k_sqr * ns_sqr; // k_{0}^{2} n_{s}^{2}
			k_sqr_ncl_sqr = k_sqr * ncl_sqr; // k_{0}^{2} n_{cl}^{2}
			k_sqr_nr_sqr = k_sqr * nr_sqr; // k_{0}^{2} n_{r}^{2}
			k_sqr_nm_sqr = k_sqr * nm_sqr; // k_{0}^{2} n_{m}^{2}

			etacs = nc_sqr / ns_sqr;
			etacr = nc_sqr / nr_sqr;
			etarcl = nr_sqr / ncl_sqr;

			// Only difference between A, B cases is search space for neff and eigenequation
			// Case A: lower = k ncl, upper = k nr, NA^{2} = nr^{2} - ncl^{2}
			// Case B: lower = k nm, upper - k nc, NA^{2} = nc^{2} - nm^{2}

			double x = nc_sqr - nm_sqr;
			double y = ns_sqr - ncl_sqr;

			na = sqrt(x); // numerical aperture

			V = (PI*d*na) / l; // V-parameter

			// predicted number of modes
			M = static_cast<int>(std::max(1.0, floor((2.0*V / PI))));

			lower = k * nm; // lower bound of search space

			upper = k * nc; // upper bound of search space

			w = k * SPEED_OF_LIGHT;

			efieldint = 0.5*EPSILON*SPEED_OF_LIGHT;

			hfieldint = 0.5*MU*SPEED_OF_LIGHT;

			// Empty the std::vector each time a new instance of the class is called
			betaE.clear();
			betaH.clear();
		}
		else {
			std::string reason = "Error: void slab_fl_neff_B::set_params(double width,double lambda,double ncore,double nsub,double nclad)\n";
			if (!c1) reason += "WG width = " + template_funcs::toString(width, 3) + " is negative\n";
			if (!c2) reason += "Wavelength = " + template_funcs::toString(lambda, 3) + " is negative\n";
			if (!c3) reason += "Cladding Index = " + template_funcs::toString(nclad, 3) + " is less than one\n";
			if (!c4) reason += "Substrate Index = " + template_funcs::toString(nsub, 3) + " is less than one\n";
			if (!c6) reason += "Rib Index = " + template_funcs::toString(nrib, 3) + " is less than one\n";
			if (!c5) reason += "Core Index = " + template_funcs::toString(ncore, 3) + " is less than cladding / substrate index\n";
			if (!c7) reason += "Core Index = " + template_funcs::toString(ncore, 3) + " is less than rib / substrate index\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_fl_neff_B::eigeneqn_b(double x, int mm, bool t)
{
	//Dispersion equation corresponding to case b from Adams
	//This means that the field oscillates in the core only	
	//This is an alternative form that produces the correct solutions

	// Case B: Field Oscillating in Core Only
	// For there to be a solution one has to have ncl < nm < nc, where nm = Max(nr,ns)

	try {

		if (k_sqr_nc_sqr > k_sqr_nm_sqr) {
			double h, p, q, r, tmp, x_sqr, v, v1, v2den, v2;

			x_sqr = template_funcs::DSQR(x);

			tmp = k_sqr_nc_sqr - x_sqr;
			h = (tmp > 0 ? sqrt(tmp) : 0.0);

			tmp = x_sqr - k_sqr_ns_sqr;
			p = (tmp > 0 ? sqrt(tmp) : 0.0);

			tmp = x_sqr - k_sqr_ncl_sqr;
			q = (tmp > 0 ? sqrt(tmp) : 0.0);
			if (!t) q *= etarcl; // multiply q by etarcl in the case of TM polarisation

			tmp = x_sqr - k_sqr_nr_sqr;
			r = (tmp > 0 ? sqrt(tmp) : 0.0);

			v = ((r - q) / (r + q)); // this includes the change to q depending on polarisation

			v1 = exp(-2.0 * r *dr); // In the limit of large dr v1 -> 0 and v2 -> 1 => Case B FLS -> TLS

			v2den = (1 + v * v1);

			v2 = (v2den > 0.0 ? (1 - v * v1) / v2den : 10.0);

			if (t) {//TE modes
				return (d * h) - atan((p / h)) - atan((r / h)*v2) - mm * PI;
			}
			else {//TM modes
				return (d * h) - atan(etacs*(p / h)) - atan(etacr * (r / h) * v2) - mm * PI;
			}
		}
		else {
			std::string reason = "Error: double slab_fl_neff_B::eigeneqn_b(double x, int mm, bool t)\n";
			reason += "Input parameters not correctly defined\n";
			reason += "k_{0}^{2} n_{c}^{2} = " + template_funcs::toString(k_sqr_nc_sqr) + ", k_{0}^{2} n_{m}^{2} = " + template_funcs::toString(k_sqr_nm_sqr) + "\n";
			return 0;
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_fl_neff_B::zbrent(double x1, double x2, double tol, bool t, int mm)
{
	//Find the root of the function between x1 and x2 using brent's method
	//The root is refined to +/- tol
	//Seemingly one of the best methods to use

	//This will be used to compute the roots of eigeneq_b
	//R. Sheehan 28 - 5 - 2010

	try {

		bool c1 = fabs(x1 - x2) > 1.0e-9 ? true : false; // cannot have x1 == x2
		bool c2 = tol > 1.0e-16 ? true : false;
		bool c3 = mm < M ? true : false;

		if (c1 && c2 && c3) {

			int iter;

			static const int ITMAX = 100;//Used in rtbis, rtflsp, rtsec, zriddr

			double a = std::min(x1, x2), b = std::max(x1, x2), c = std::max(x1, x2), d, e, min1, min2;
			double fc, p, q, r, s, tol1, xm;
			double fa = eigeneqn_b(a, mm, t), fb = eigeneqn_b(b, mm, t);

			if ((fa>0.0 && fb>0.0) || (fa<0.0 && fb<0.0)) {
				std::cerr << "Root must be bracketed in zbrent\n";
			}
			fc = fb;
			for (iter = 1; iter <= ITMAX; iter++) {
				if ((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0)) {
					c = a;
					fc = fa;
					e = d = b - a;
				}
				if (fabs(fc)<fabs(fb)) {
					a = b;
					b = c;
					c = a;
					fa = fb;
					fb = fc;
					fc = fa;
				}
				tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
				xm = 0.5*(c - b);
				if (fabs(xm) <= tol1 || fb == 0.0) return b;
				if(fabs(xm)<=tol1 || fb==0.0){
					std::cout<<"Brent's Method converged in "<<iter<<" iterations\n";
					return b;
				}
				if (fabs(e) >= tol1 && fabs(fa)>fabs(fb)) {
					s = fb / fa;
					if (a == c) {
						p = 2.0*xm*s;
						q = 1.0 - s;
					}
					else {
						q = fa / fc;
						r = fb / fc;
						p = s * (2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
						q = (q - 1.0)*(r - 1.0)*(s - 1.0);
					}
					if (p>0.0) q = -q;
					p = fabs(p);
					min1 = 3.0*xm*q - fabs(tol1*q);
					min2 = fabs(e*q);
					if (2.0*p<std::min(min1, min2)) {
						e = d;
						d = p / q;
					}
					else {
						d = xm;
						e = d;
					}
				}
				else {
					d = xm;
					e = d;
				}
				a = b;
				fa = fb;
				if (fabs(d)>tol1) {
					b += d;
				}
				else {
					b += template_funcs::SIGN(tol1, xm);
				}
				fb = eigeneqn_b(b, mm, t);
			}
			std::cerr << "Maximum number of iterations exceeded in zbrent\n";
			return 0.0;
		}
		else {
			std::string reason = "Error: double slab_fl_neff_B::zbrent(double x1,double x2,double tol,bool t,int mm)\n";
			if (!c1) reason += "Cannot have x1 = x2\nx1 = " + template_funcs::toString(x1) + ", x2 = " + template_funcs::toString(x1) + "\n";
			if (!c2) reason += "Desired tolerance is less than smallest allowable EPS\n";
			if (!c3) reason += "mm >= M\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_fl_neff_B::neff_search(bool mode)
{
	// Compute the waveguide mode effective indices based on given polarisation and wg_type

	// Case B: Field Oscillating in Core Only
	// For there to be a solution one has to have ncl < nm < nc, where nm = Max(nr,ns)

	// Case B: lower = k nm, upper - k nc, NA^{2} = nc^{2} - nm^{2}

	try {

		if (lower < upper) {

			int m;
			double b;

			std::vector<double> vec;

			for (m = 0; m < M; m++) {
				b = zbrent(lower, upper, EPS, mode, m);

				if (b>lower && b<upper) {
					vec.push_back(b);
				}

			}

			if (mode) {
				betaE = vec;
			}
			else {
				betaH = vec;
			}

			vec.clear();

		}
		else {
			std::string reason = "Error: void slab_fl_neff_B::neff_search(bool mode, bool wg_type)\n";
			reason += "Search range is not correctly defined\n";
			reason += "lower = " + template_funcs::toString(lower, 4) + ", upper = " + template_funcs::toString(upper, 4) + "\n";
			throw std::range_error(reason);
		}
	}
	catch (std::range_error &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_fl_neff_B::output_disp_eqn_b(bool mode, std::string& storage_directory)
{
	// Output the case A four-layer slab dispersion equation
	// R. Sheehan 18 - 12 - 2020

	try {
		if (lower < upper) {
			int N = 101;

			// create array to hold dispersion equations for all the modes
			std::vector< std::vector< double > > mat;

			mat.resize(M + 2); // We want to output M solutions plus the corresponding positions

			// store the positions at which the eigeneqn_a will be evaluated
			double dx = (upper - lower) / (static_cast<double>(N - 1));

			double xi = lower;

			mat[0].resize(N + 2);

			for (int j = 1; j <= N; j++) {
				mat[0][j] = xi;
				xi += dx;
			}

			// compute the values of the dispersion equation for each of the expected modes
			for (int i = 1; i <= M; i++) {
				mat[i].resize(N + 2);
				for (int j = 1; j <= N; j++) {
					mat[i][j] = eigeneqn_b(mat[0][j], i - 1, mode);
				}
			}

			// Output the computed dispersion equations
			std::string pol = (mode ? "TE" : "TM");
			std::string filename = storage_directory + pol + "_Disp_Eqns_FL_B.txt";
			std::ofstream write;

			write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);

			if (write.is_open()) {

				for (int i = 1; i <= N; i++) {
					for (int j = 0; j < (M + 1); j++)
						if (j == M)
							write << std::setprecision(20) << mat[j][i];
						else
							write << std::setprecision(20) << mat[j][i] << " , ";
					write << "\n";
				}

				write.close();
			}

			mat.clear();
		}
		else {
			std::string reason = "Error: void slab_fl_neff_B::output_disp_eqn(bool mode)\n";
			reason += "Search range is not correctly defined\n";
			reason += "lower = " + template_funcs::toString(lower, 4) + ", upper = " + template_funcs::toString(upper, 4) + "\n";
			throw std::range_error(reason);
		}
	}
	catch (std::range_error& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_fl_neff_B::report(bool mode)
{
	// print the computed results to screen
	std::cout << "Output for the slab waveguide calculator\n";
	std::cout << "\n";
	std::cout << "Waveguide width = " << d << " microns\n";
	std::cout << "Rib width = " << dr << " microns\n";
	std::cout << "Wavelength = " << l << " microns\n";
	std::cout << "Wavenumber = " << k << " inverse microns\n";
	std::cout << "\n";
	std::cout << "Core Index = " << nc << "\n";
	std::cout << "Substrate Index = " << ns << "\n";
	std::cout << "Rib Index = " << nr << "\n";
	std::cout << "Cladding Index = " << ncl << "\n";
	std::cout << "\n";
	std::cout << "Numerical Aperture = " << na << "\n";
	std::cout << "V-Parameter = " << V << "\n";
	std::cout << "Number of Modes = " << M << "\n";
	std::cout << "\n";

	std::string pol = (mode ? "TE" : "TM");

	std::cout << pol << " Modes\n";
	std::cout << "There are " << nbeta(mode) << " calculated " + pol + " modes\n";

	for (int i = 0; i < nbeta(mode); i++) {
		std::cout << pol + "_{" << i + 1 << "} = " << std::setprecision(10) << _beta(i, mode) << " , n_eff_{" << i + 1 << "} = " << std::setprecision(10) << (_beta(i, mode) / k) << "\n";
	}
	std::cout << "\n";
}

int slab_fl_neff_B::get_nmodes(bool mode)
{
	// return the number of computed propagation constants for a given polarisation
	// it can be the case that the predicted value M is greater than the actual number of modes in a waveguide

	return nbeta(mode);
}

double slab_fl_neff_B::get_neff(int i, bool mode)
{
	// return the i^{th} effective index for a given polarisation
	// by convention mode == true => TE modes, mode == false => TM modes

	try {
		if (nbeta(mode) > 0) {
			if (i<0 || i > nbeta(mode)) {
				throw std::range_error("Error: double slab_fl_neff_B::get_neff(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (mode) {
					return betaE[i] / k;
				}
				else {
					return betaH[i] / k;
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_fl_neff_B::get_neff(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error& e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// slab_fl_mode_B is used to compute the mode profile in the four layer slab case B

slab_fl_mode_B::slab_fl_mode_B()
{
	// Default
}

slab_fl_mode_B::slab_fl_mode_B(double width, double rib_width, double lambda, double ncore, double nsub, double nclad, double nrib)
{
	set_params(width, rib_width, lambda, ncore, nsub, nclad, nrib);
}

void slab_fl_mode_B::compute_neff(bool mode)
{
	// Compute the effective indices of the slab waveguide
	neff_search(mode);

	std::string pol = (mode ? "TE" : "TM");

	std::cout << "There are " << nbeta(mode) << " calculated " + pol + " modes\n";
	for (int i = 0; i<static_cast<int>(nbeta(mode)); i++) {
		std::cout << "beta[" << i + 1 << "] = " << std::setprecision(9) << _beta(i, mode) << ", "<< _beta(i, mode) / k << "\n";
	}
	std::cout << "\n";
}

double slab_fl_mode_B::phase(int i, bool t)
{
	// Phase of mode in type B four layer slab waveguide

	// phase term may not be valid for this case

	try {
		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_fl_mode_B::phase(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					double qq = q(i, t);
					double rr = rB(i, t); 
					double arg = qq / rr; 
					double atval = fabs(arg) > 1 ? PI_2 : atanh(arg); 
					return -1.0*( ( dr * rr ) + atval );
				}
				else {
					double qq = q(i, t);
					double rr = rB(i, t);
					double arg = ( (etarcl * qq) / rr);
					double atval = fabs(arg) > 1 ? PI_2 : atanh(arg);
					return -1.0*( ( dr * rr ) + atval);
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_fl_mode_B::phase(int i, bool t)s\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_fl_mode_B::TE_TM(double x, int i, bool mode) 
{
	// shape of waveguide mode
	// Function which defines the shape of the modes for Type B four layer slab
	// This method uses the stored computed propagation constants
	// Assumes that the slab wg has core region defined on -W < x < 0 and 0 < x < D

	// https://www.smbc-comics.com/comic/cure

	try {		
		if (nbeta(mode) > 0) {
			if (i<0 || i > nbeta(mode)) {
				throw std::range_error("Error: double slab_fl_mode_B::TE_TM(double x,int i,bool mode)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				double hh = h(i, mode); 
				double rr = rB(i, mode); 
				double wh = d * hh;
				double qq = q(i, mode);
				double pp = p(i, mode); 
				double ph = phase(i, mode); 
				double rrhh = mode ? (rr / hh) * tanh(ph) : (rr / hh) * tanh(ph) * etacr;
				if (x < -d) {
					// x < -W
					return ( ( cos(wh) - ( rrhh * sin(wh) ) )*exp( pp * ( x + d ) ) );
				}
				else if (x >= -d && x <= 0) {
					// -W < x < 0
					return ( cos( hh * x ) + rrhh * sin( hh * x ) ); 
				}
				else if (x > 0 && x <= dr) {
					// 0 < x < D
					return ( cosh( (rr * x) + ph) / cosh(ph) ); 
				}
				else {
					// x > D
					return ( ( cosh( (rr * dr) + ph) / cosh(ph) )*exp( qq * ( dr - x ) ) ); 
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_fl_mode_B::TE_TM(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_fl_mode_B::output_modes(bool mode, int N, double Lx, std::string& storage_directory)
{
	// This function will calculate the solutions corresponding to each mode and output them to a file for a single polarisation

	double xi;

	std::string filename;

	std::string pol = (mode ? "TE" : "TM");

	std::ofstream write;

	if (nbeta(mode) > 0) {

		std::vector< std::vector< double > > mat;

		mat.resize(nbeta(mode) + 2); // We want to output M solutions plus the corresponding positions

		double dx = Lx / (static_cast<double>(N - 1));

		xi = -0.5 * Lx;

		mat[0].resize(N + 2);

		for (int j = 1; j <= N; j++) {
			mat[0][j] = xi;
			xi += dx;
		}

		for (int i = 1; i <= nbeta(mode); i++) {
			mat[i].resize(N + 2);
			for (int j = 1; j <= N; j++) {
				mat[i][j] = TE_TM(mat[0][j], i - 1, mode);
			}
		}

		// Output all the modes to a single file
		filename = storage_directory + pol + "_Mode_Profiles_FL_B.txt";
		write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);

		if (write.is_open()) {

			for (int i = 1; i <= N; i++) {
				for (int j = 0; j < (nbeta(mode) + 1); j++)
					if (j == nbeta(mode))
						write << std::setprecision(20) << mat[j][i];
					else
						write << std::setprecision(20) << mat[j][i] << " , ";
				write << "\n";
			}

			write.close();

		}

		output_disp_eqn_b(mode, storage_directory);
	}
}

// definitions for the coupled slab waveguide class

coupled_slab_tl_neff::coupled_slab_tl_neff()
{
	// Default
	slab_sep = coupling_coeff = L_coupling = 0.0; 
}

coupled_slab_tl_neff::coupled_slab_tl_neff(double separation, double width, double lambda, double ncore, double nsub)
{
	// Assign value to the slab parameters

	set_params(separation, width, lambda, ncore, nsub); 
}

void coupled_slab_tl_neff::set_params(double separation, double width, double lambda, double ncore, double nsub)
{
	// Assign values to the parameters for the coupled slab structure

	try {	
		if (separation > 0.0) {
			slab_tl_neff::set_params(width, lambda, ncore, nsub, nsub);

			slab_sep = separation;
		}
		else {
			throw std::invalid_argument("Error: coupled_slab_tl_neff::set_params(double separation, double width, double lambda, double ncore, double nsub)\nSeparation < 2 Width\n");
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	} 
}

double coupled_slab_tl_neff::compute_coupling_coeff(bool mode)
{
	// Compute the coupling coefficient for the coupled waveguides
	// theory taken from Okamoto, "Fundamentals of Optical Waveguides"
	// returned value is in units of um^{-1}
	// R. Sheehan 1 - 10 - 2018

	try {
		neff_search(mode);

		if (nbeta(mode) > 0) {

			double tt = d / 2; 
			double hh = h(0, mode); 
			double pp = p(0, mode);
			double vs = template_funcs::DSQR(k * na);
			double num = template_funcs::DSQR(hh * pp); 
			double denom = _beta(0, mode)*(1.0 + ( tt * pp ) )*vs;
			double arg = -1.0 * pp * (slab_sep - d); 
			
			if (denom > 0) {
				return ( (num / denom)* exp(arg) );
			}
			else {
				throw std::runtime_error("Error: double coupled_slab_tl_neff::compute_coupling_coeff(bool mode)\nDivision by zero occurred\n");
				return 0.0; 
			}
		}
		else {
			throw std::runtime_error("Error: double coupled_slab_tl_neff::compute_coupling_coeff(bool mode)\nNo modes computed for this waveguide\n"); 
			return 0;
		}
	}
	catch (std::runtime_error &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// definitions for the coupled slabs class
coupled_slabs::coupled_slabs()
{
	// default constructor
	wg_defined = false; coeffs_defined = false; 
	pol = TE; // going to assume TE polarisation for all calcs for simplicity
	msize = 2; 
	min_pitch = de_A = de_B = omega = wavenum = cw1 = cw2 = WA = WB = n_core_A = n_core_B = n_sub = wavel = 0.0;
	CAA = CBB = CAB = norm = KAA = KBB = KAB = KBA = ga = gb = kab = kba = async = 0.0; 
	bp = bm = phi = psi = del = Lc = 0.0; 
}

coupled_slabs::coupled_slabs(double W1, double W2, double lambda, double ncore1, double ncore2, double nsub)
{
	// primary constructor
	set_params(W1, W2, lambda, ncore1, ncore2, nsub);
}

void coupled_slabs::set_params(double W1, double W2, double lambda, double ncore1, double ncore2, double nsub)
{
	// Method for assigning values to each of the WG
	
	try {
		bool c1 = W1 > 0.0 ? true : false; 
		bool c2 = W2 > 0.0 ? true : false; 
		bool c4 = lambda > 0.0 ? true : false; 
		bool c5 = nsub > 0.0 ? true : false; 
		bool c6 = ncore1 > nsub ? true : false; 
		bool c7 = ncore2 > nsub ? true : false; 
		bool c10 = c1 && c2 && c4 && c5 && c6 && c7; 
		
		if (c10) {

			// store the parameters locally for an easy life
			WA = W1; WB = W2; n_core_A = ncore1; n_core_B = ncore2; n_sub = nsub; wavel = lambda; 

			min_pitch = 0.5 * (WA + WB);
			
			// w = k * SPEED_OF_LIGHT; 
			wavenum = Two_PI / wavel; 

			omega = ( SPEED_OF_LIGHT * wavenum ); 
			
			cw1 = ( EPSILON * omega ) / 4.0; // constant required in calculation
			
			cw2 = 2.0 * omega * MU; // constant required in calculation
			
			//de_A = ( cw1 * ( template_funcs::DSQR(n_core_B) - template_funcs::DSQR(n_sub) ) ) ;
			de_A = ( cw1 * template_funcs::SQR_DIFF(n_core_B, n_sub) ) ;

			//de_B = ( cw1 * ( template_funcs::DSQR(n_core_A) - template_funcs::DSQR(n_sub) ) ) ;
			de_B = ( cw1 * template_funcs::SQR_DIFF(n_core_A, n_sub) ) ;

			WGA.set_params(W1, lambda, ncore1, nsub, nsub); // assign the parameters to the slab object
			
			WGA.compute_neff(pol); // compute the effective indices of that slab
			
			WGB.set_params(W2, lambda, ncore2, nsub, nsub); // assign the parameters to the slab object
			
			WGB.compute_neff(pol); // compute the effective indices of that slab

			wg_defined = true; 
		}
		else {
			std::string reason = "Error: void coupled_slabs::set_params(double W1, double W2, double separ, double lambda, double ncore1, double ncore2, double nsub)\n";
			if (!c1) reason += "WGA width = " + template_funcs::toString(W1, 3) + " is negative\n";
			if (!c2) reason += "WGB width = " + template_funcs::toString(W2, 3) + " is negative\n";
			if (!c4) reason += "Wavelength = " + template_funcs::toString(lambda, 3) + " is negative\n";
			
			if (!c5) reason += "Substrate Index = " + template_funcs::toString(nsub, 3) + " is too small\n";
			if (!c6) reason += "WGA Index = " + template_funcs::toString(ncore1, 3) + " is too small\n";			
			if (!c7) reason += "WGB Index = " + template_funcs::toString(ncore2, 3) + " is too small\n";			

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument & e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void coupled_slabs::compute_coefficients(double pitch, bool loud)
{
	// Compute the various coupling coefficients

	try {
		
		// Start by checking that you can plot the mode profiles together with the correct separations
		// Then check that you can integrate the mode profiles over an appropriate range
		// Then test calculation of C and K

		bool c3 = pitch > min_pitch ? true : false;

		if (c3 && wg_defined) {

			bool loud_integral = false; 
			bool scale_integral = true; 

			// overlap integrals
			CAA = integrate_CAA(pitch, scale_integral, loud_integral);
			CBB = integrate_CBB(pitch, scale_integral, loud_integral);
			CAB = integrate_CAB(pitch, scale_integral, loud_integral);
			norm = CAA + CBB; 
			CAB /= norm; 
			
			// coupling coefficients
			KAA = integrate_KAA(pitch, scale_integral, loud_integral);
			KBB = integrate_KBB(pitch, scale_integral, loud_integral);
			KAB = integrate_KAB(pitch, scale_integral, loud_integral);
			KBA = integrate_KBA(pitch, scale_integral, loud_integral);

			KAA /= norm; KBB /= norm; KAB /= norm; KBA /= norm; 

			// propagation matrix elements
			double bA, bB, numer; 
			numer = 1.0 - template_funcs::DSQR(CAB); 
			bA = WGA.get_neff(0, pol) * wavenum; 
			bB = WGB.get_neff(0, pol) * wavenum; 			
			
			ga = bA + ((KAA - CAB * KBA) / numer); 
			gb = bB + ((KBB - CAB * KAB) / numer); 
			kab = ((KAB - CAB * KBB) / numer); 
			kba = ((KBA - CAB * KAA) / numer); 

			// asynchronism factor
			async = (gb - ga) / (2.0 * sqrt(kab * kba)); 

			// propagation constants in coupled waveguides
			phi = 0.5 * (ga + gb); 
			del = 0.5 * (gb - ga); 
			psi = sqrt(template_funcs::DSQR(del) + kab * kba); 
			bp = phi + psi; 
			bm = phi - psi; 
			Lc = PI / fabs(bm - bp); 

			double t1 = del + psi; 
			double t2 = psi - del; 

			// compute the elements of the matrices V and Vinv
			V = vecut::zero_cmat(msize, msize); 
			Vinv = vecut::zero_cmat(msize, msize);

			V[0][0] = kab; V[0][1] = kab; 
			V[1][0] = -1.0 * t1; V[1][1] = t2; 

			double detV = 2.0 * kab * psi; 
			Vinv[0][0] = t2 / detV; Vinv[1][1] = kab / detV; 
			Vinv[1][0] = t1 / detV; Vinv[0][1] = -1.0 * Vinv[1][1]; 

			coeffs_defined = true; 

			if (loud) {
				std::cout << "Computed Coefficients\n";
				std::cout << "Overlap Integrals\n";
				std::cout << "C_{aa} = " << CAA << "\n";
				std::cout << "C_{bb} = " << CBB << "\n";
				std::cout << "C_{ab} = " << CAB << "\n\n";
				//std::cout << "C / (C_{aa} + C_{bb}) = " << CAB / norm << "\n\n";

				std::cout << "Coupling Coefficients\n";
				std::cout << "K_{aa} = " << KAA << "\n";
				std::cout << "K_{bb} = " << KBB << "\n";
				std::cout << "K_{ab} = " << KAB << "\n";
				std::cout << "K_{ba} = " << KBA << "\n\n";

				std::cout << "Overlap-integral-coupling-coefficient relation\n";
				std::cout << "(beta_{b} - beta_{a}) C = " << (bB - bA) * CAB << "\n";
				std::cout << "K_{ba} - K_{ab} = " << KBA - KAB << "\n\n";

				std::cout << "Propagation Matrix Elements\n";
				std::cout << "g_{a} = " << ga << "\n";
				std::cout << "g_{b} = " << gb << "\n";
				std::cout << "k_{ab} = " << kab << "\n";
				std::cout << "k_{ba} = " << kba << "\n\n";

				std::cout << "Propagation Matrix Element Relation\n";
				std::cout << "k_{ab} - k_{ba} = " << kab - kba << "\n";
				std::cout << "(g_{a} - g_{b}) C = " << (ga - gb) * CAB << "\n";
				std::cout << "| (k_{ab} - k_{ba}) - (g_{a} - g_{b}) C | = " << fabs(kab - kba) - fabs((ga - gb) * CAB) << "\n\n";

				std::cout << "Asynchronism\n";
				std::cout << "d = " << async << "\n\n";

				std::cout << "Propagation constants\n"; 
				std::cout << "bp = " << bp << " um^{-1}\n"; 
				std::cout << "bm = " << bm << "um^{-1}\n\n"; 

				std::cout << "Coupling Length\n"; 
				std::cout << "Lc = " << Lc << "um\n\n"; 

				std::cout << "Matrix V\n"; 
				vecut::print_cmat(V); 
				std::cout << "\n"; 

				std::cout << "det(V) = " << detV << "\n\n"; 

				std::cout << "Matrix V^{-1}\n";
				vecut::print_cmat(Vinv);
				std::cout << "\n";
			}
		}
		else {
			std::string reason = "Error: void coupled_slabs::compute_coefficients(double pitch)\n";
			if (!wg_defined) reason += "WG parameters are not defined\n"; 
			if (!c3) reason += "WG pitch = " + template_funcs::toString(pitch, 3) + " is too small\n";
			if (!c3) reason += "min pitch = " + template_funcs::toString(min_pitch, 3) + "\n";
			
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument & e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void coupled_slabs::output_modes(double pitch)
{
	// output the modes in the same file separated by their pitch
	
	try {

		bool c3 = pitch > min_pitch ? true : false;

		if (c3 && wg_defined) {
			int N = 501; 
			double Lx = 2.5 * (WA + pitch + WB); // total length of simulation region
			double xmid = 0.5 * pitch + 0.25 * (WA + WB); // midpoint of simulation region
			double x0 = xmid - 0.5 * Lx; // start computing solutions here
			double dx = Lx / (N - 1); 
			double vA, vB; 

			std::string filename = "Coupled_Mode_Profiles.txt";
			std::ofstream write;

			write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);

			if (write.is_open()) {
				
				for (int i = 0; i < N; i++) {
					vA = WGA.TE_TM(x0, 0, pol); 
					vB = WGB.TE_TM(x0 - pitch, 0, pol); 
					write << std::setprecision(10) << x0 << " , " << vA << " , " << vB << "\n"; 

					// store the field values for later use
					Afield.push_back(std::complex<double>(vA, 0.0)); 
					Bfield.push_back(std::complex<double>(vB, 0.0)); 
					pos.push_back(x0); 

					x0 += dx; 
				}

				write.close(); 
			}
		}
		else {
			std::string reason = "Error: void coupled_slabs::output_modes(double pitch)\n";
			if (!wg_defined) reason += "WG parameters are not defined\n";
			if (!c3) reason += "WG pitch = " + template_funcs::toString(pitch, 3) + " is too small\n";
			if (!c3) reason += "min pitch = " + template_funcs::toString(min_pitch, 3) + "\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

//double coupled_slabs::integrate_modes(int integrand, double x1, double x2, double pitch)
//{
//	// dynamic binding version is much simpler to implement and understand and debug
//	// R. Sheehan 9 - 9 - 2020
//
//	// compute the integral of the modes based on the choice of integrand
//	// integrand == 0 => pq == aa => E_{ay} E_{ay}
//	// integrand == 1 => pq == bb => E_{by} E_{by}
//	// integrand == 2 => pq == ab => E_{ay} E_{by}
//	// integrand == 3 => pq == ba => E_{by} E_{ay}
//
//	try {
//		bool c1 = fabs(x2 - x1) > EPS ? true : false;
//		bool c2 = integrand > -1 && integrand < 4 ? true : false; 
//		bool c3 = pitch > 0.5 * (WA + WB) ? true : false;
//		bool c10 = c1 && c3; 
//
//		if (c10 && wg_defined) {
//			int N = 501;
//			double Lx = 2.5 * (WA + pitch + WB); // total length of simulation region
//			double xmid = 0.5 * pitch + 0.25 * (WA + WB); // midpoint of simulation region
//			double x0 = xmid - 0.5 * Lx; // start computing solutions here
//			double xoff = x0 - pitch; // use the offset value to evaluate mode solution in WGB
//			double dx = Lx / (N - 1);
//			double first, last, term, integral, sum = 0.0; 
//
//			// there's no point in actually having two function evaluations is there?
//			// use a switch to evaluate both mode solution values only when required
//			// otherwise compute one value and square it
//			bool evaluate_each_mode = false; 
//
//			// Use a pointer to the slab waveguide object to access the data stored therein
//			// this will obviate the need to repeatedly decide which integrand must be specified
//			// during the loop for computing the integrals
//
//			slab_tl_mode* md1 = NULL; 
//			slab_tl_mode* md2 = NULL;
//			
//			if (integrand == 0) {
//				// integrand == 0 => pq == aa => E_{ay} E_{ay}
//				md1 = &WGA; 
//				md2 = &WGA; 
//				evaluate_each_mode = false; // tell alg. to evaluate value of field in WGA only and square it
//			}
//			else if (integrand == 1) {
//				// integrand == 1 => pq == bb => E_{by} E_{by}
//				md1 = &WGB;
//				md2 = &WGB;
//				evaluate_each_mode = false; // tell alg. to evaluate value of field in WGB only and square it
//			}
//			else if (integrand == 2 || integrand == 3) {
//				// integrand == 2 => pq == ab => E_{ay} E_{by}
//				// integrand == 3 => pq == ba => E_{by} E_{ay}
//				md1 = &WGA;
//				md2 = &WGB;
//				evaluate_each_mode = true; // tell alg. to evaluate value of field in both WG
//			}
//			else {
//				// you done fucked up bro!
//				// you should not be here!
//				// program will have already crashed 
//				// if integrand < 0 || integrand > 3
//			}
//
//			// you still end up having to make decisions at each step of the loop 
//			// the question is whether or not it is more efficient to 
//			// a) have a single function requiring multiple decisions at each step
//			// b) have multiple functions that essentially do the same thing but each function operates on a different integrand
//			// Is dynamic binding more efficient than if-else? 
//			// The answer seems to be that from an efficiency POV it doesn't matter, switches and dynamic binding are equally fast
//			// However, from code design and maintenance and scalability POV then dynamic binding is definitely better
//			// because it becomes much easier to add cases as the need arises
//			// https://stackoverflow.com/questions/2681337/dynamical-binding-or-switch-case
//			// I would agree that it will be easier to debug the individual cases so I might go that way. 
//
//			// Evaluate the integral using the trapezoidal rule
//
//			// compute first term in the sum
//			if (evaluate_each_mode) {
//				// you must be computing an integral of E_{ay} E_{by}
//				first = md1->TE_TM(x0, 0, pol) * md2->TE_TM(xoff, 0, pol);
//			}
//			else {
//				// you must be computing an integral of E_{ay} E_{ay} or E_{by} E_{by} 
//				first = template_funcs::DSQR(integrand == 0 ? md1->TE_TM(x0, 0, pol) : md1->TE_TM(xoff, 0, pol));
//			}
//
//			// update position
//			x0 += dx;
//			xoff += dx;
//
//			// sum over the middle terms 
//			for (int i = 1; i < N-1; i++) {
//				
//				if (evaluate_each_mode) {
//					// you must be computing an integral of E_{ay} E_{by}
//					term = md1->TE_TM(x0, 0, pol) * md2->TE_TM(xoff, 0, pol);
//				}
//				else {
//					// you must be computing an integral of E_{ay} E_{ay} or E_{by} E_{by} 
//					term = template_funcs::DSQR( integrand == 0 ? md1->TE_TM(x0, 0, pol) : md1->TE_TM(xoff, 0, pol) ); 
//				}
//				
//				sum += term;
//
//				// update position
//				x0 += dx;
//				xoff += dx; 
//			}
//
//			// compute last term in the sum
//			if (evaluate_each_mode) {
//				// you must be computing an integral of E_{ay} E_{by}
//				last = md1->TE_TM(x0, 0, pol) * md2->TE_TM(xoff, 0, pol);
//			}
//			else {
//				// you must be computing an integral of E_{ay} E_{ay} or E_{by} E_{by} 
//				last = template_funcs::DSQR( integrand == 0 ? md1->TE_TM(x0, 0, pol) : md1->TE_TM(xoff, 0, pol));
//			}
//
//			integral = 0.5 * dx * (first + last + 2.0 * sum); 
//
//			return integral; 
//		}
//		else {
//			return 0; 
//
//			std::string reason = "Error: double coupled_slabs::integrate_modes(double x1, double x2, double pitch)\n";
//			if (!wg_defined) reason += "WG parameters are not defined\n";
//			if (!c1) reason += "Integration Limits are not correct\n";
//			if (!c2) reason += "Integrand not correctly specified\n";
//			if (!c3) reason += "WG pitch = " + template_funcs::toString(pitch, 3) + " is too small\n";
//
//			throw std::invalid_argument(reason);
//		}
//	}
//	catch (std::invalid_argument& e) {
//		useful_funcs::exit_failure_output(e.what());
//		exit(EXIT_FAILURE);
//	}
//}

double coupled_slabs::integrate_CAA(double pitch, bool scale, bool loud)
{
	// compute the overlap integral of E_{ay} E_{ay}
	
	try {
		//bool c1 = fabs(x2 - x1) > EPS ? true : false;
		//bool c3 = pitch > 0.5 * (WA + WB) ? true : false;
		//bool c10 = c1 && c3;

		bool c3 = true; 

		if (c3 && wg_defined) {
			int N = 501;
			double Lx = 2.5 * (WA + pitch + WB); // total length of simulation region
			double xmid = 0.5 * pitch + 0.25 * (WA + WB); // midpoint of simulation region
			double x0 = xmid - 0.5 * Lx; // start computing solutions here
			double dx = Lx / (N - 1);
			double first, last, term, integral, sum = 0.0;

			std::string filename = "integrate_CAA_field_values.txt";
			std::ofstream write;

			if (loud) {
				write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);
			}

			// Evaluate the integral using the trapezoidal rule

			// compute first term in the sum
			first = template_funcs::DSQR( WGA.TE_TM(x0, 0, pol) );

			// update position
			x0 += dx;

			// sum over the middle terms 
			for (int i = 1; i < N - 1; i++) {

				if (loud) {
					write << std::setprecision(10) << x0 << " , " << WGA.TE_TM(x0, 0, pol) << " , " << WGB.TE_TM(x0 , 0, pol) << "\n";
				}

				term = template_funcs::DSQR( WGA.TE_TM(x0, 0, pol) );

				sum += term;

				// update position
				x0 += dx;
			}

			if (loud) {
				write.close(); 
			}

			// compute last term in the sum
			last = template_funcs::DSQR( WGA.TE_TM(x0, 0, pol) );

			integral = 0.5 * dx * (first + last + 2.0 * sum);

			if(scale) integral *= ( WGA.get_neff(0, pol) * wavenum ) / cw2; 

			return integral;
		}
		else {
			return 0;

			std::string reason = "Error: double coupled_slabs::integrate_CAA(double pitch)\n";
			if (!wg_defined) reason += "WG parameters are not defined\n";
			//if (!c1) reason += "Integration Limits are not correct\n";
			if (!c3) reason += "WG pitch = " + template_funcs::toString(pitch, 3) + " is too small\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double coupled_slabs::integrate_CBB(double pitch, bool scale, bool loud)
{
	// compute the overlap integral of E_{by} E_{by}

	try {
		//bool c1 = fabs(x2 - x1) > EPS ? true : false;
		//bool c3 = pitch > 0.5 * (WA + WB) ? true : false;
		//bool c10 = c1 && c3;

		bool c3 = true; 

		if (c3 && wg_defined) {
			int N = 501;
			double Lx = 2.5 * (WA + pitch + WB); // total length of simulation region
			double xmid = 0.5 * pitch + 0.25 * (WA + WB); // midpoint of simulation region
			double x0 = xmid - 0.5 * Lx; // start computing solutions here
			double xoff = x0 - pitch; // use the offset value to evaluate mode solution in WGB
			double dx = Lx / (N - 1);
			double first, last, term, integral, sum = 0.0;

			std::string filename = "integrate_CBB_field_values.txt";
			std::ofstream write;

			if (loud) {
				write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);
			}

			// Evaluate the integral using the trapezoidal rule

			// compute first term in the sum
			first = template_funcs::DSQR(WGB.TE_TM(xoff, 0, pol));

			// update position
			xoff += dx;

			// sum over the middle terms 
			for (int i = 1; i < N - 1; i++) {

				if (loud) {
					write << std::setprecision(10) << xoff << " , " << WGB.TE_TM(xoff, 0, pol) << " , " << WGB.TE_TM(xoff, 0, pol) << "\n";
				}

				term = template_funcs::DSQR(WGB.TE_TM(xoff, 0, pol));

				sum += term;

				// update position
				xoff += dx;
			}

			if (loud) {
				write.close();
			}

			// compute last term in the sum
			last = template_funcs::DSQR(WGB.TE_TM(xoff, 0, pol));

			integral = 0.5 * dx * (first + last + 2.0 * sum);

			if(scale) integral *= (WGB.get_neff(0, pol) * wavenum) / cw2;

			return integral;
		}
		else {
			return 0;

			std::string reason = "Error: double coupled_slabs::integrate_CBB(double pitch)\n";
			if (!wg_defined) reason += "WG parameters are not defined\n";
			//if (!c1) reason += "Integration Limits are not correct\n";
			if (!c3) reason += "WG pitch = " + template_funcs::toString(pitch, 3) + " is too small\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double coupled_slabs::integrate_CAB(double pitch, bool scale, bool loud)
{
	// compute the overlap integral of E_{ay} E_{by}
	// This corresponds to constant C in the notes

	try {
		//bool c3 = pitch > 0.5 * (WA + WB) ? true : false;

		bool c3 = true; 

		if (c3 && wg_defined) {
			int N = 501;
			double Lx = 2.5 * (WA + pitch + WB); // total length of simulation region
			double xmid = 0.5 * pitch + 0.25 * (WA + WB); // midpoint of simulation region
			double x0 = xmid - 0.5 * Lx; // start computing solutions here
			double xoff = x0 - pitch; // use the offset value to evaluate mode solution in WGB
			double dx = Lx / (N - 1);
			double first, last, term, integral, sum = 0.0;

			std::string filename = "integrate_CAB_field_values.txt";
			std::ofstream write;

			if (loud) {
				write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);
			}

			// Evaluate the integral using the trapezoidal rule

			// compute first term in the sum
			first = ( WGA.TE_TM(x0, 0, pol) * WGB.TE_TM(xoff, 0, pol) );

			// update position
			x0 += dx;
			xoff += dx;

			// sum over the middle terms 
			for (int i = 1; i < N - 1; i++) {

				if (loud) {
					write << std::setprecision(10) << x0 << " , " << WGA.TE_TM(x0, 0, pol) << " , " << WGB.TE_TM(xoff, 0, pol) << "\n";
				}

				term = (WGA.TE_TM(x0, 0, pol) * WGB.TE_TM(xoff, 0, pol));

				sum += term;

				// update position
				x0 += dx;
				xoff += dx;
			}

			if (loud) {
				write.close();
			}

			// compute last term in the sum
			last = (WGA.TE_TM(x0, 0, pol) * WGB.TE_TM(xoff, 0, pol));

			integral = 0.5 * dx * (first + last + 2.0 * sum); 

			if(scale) integral *= ( 0.5 * ( WGA.get_neff(0, pol) + WGB.get_neff(0, pol) ) * wavenum) / cw2; 

			integral *= 2.0; 

			return integral;
		}
		else {
			return 0;

			std::string reason = "Error: double coupled_slabs::integrate_CAB(double pitch)\n";
			if (!wg_defined) reason += "WG parameters are not defined\n";
			//if (!c1) reason += "Integration Limits are not correct\n";
			if (!c3) reason += "WG pitch = " + template_funcs::toString(pitch, 3) + " is too small\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double coupled_slabs::integrate_KAA(double pitch, bool scale, bool loud)
{
	// compute the coupling integral of E_{ay} E_{ay}
	// range of integration is D-b <= x <= D+b

	try {
		bool c3 = pitch > 0.5 * (WA + WB) ? true : false;

		if (c3 && wg_defined) {
			int N = 501;
			double Lx = 2.5 * (WA + pitch + WB); // total length of simulation region
			double dx = Lx / (N - 1);
			double x0 = pitch - 0.5 * WB; 
			double xend = pitch + 0.5 * WB; 
			double first, last, term, integral, sum = 0.0;

			std::string filename = "integrate_KAA_field_values.txt";
			std::ofstream write;

			if (loud) {
				write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);
			}

			// Evaluate the integral using the trapezoidal rule

			// compute first term in the sum
			first = template_funcs::DSQR(WGA.TE_TM(x0, 0, pol));

			// update position
			x0 += dx;

			// sum over the middle terms 
			while(x0 < xend) {

				if (loud) {
					write << std::setprecision(10) << x0 << " , " << WGA.TE_TM(x0, 0, pol) << " , " << WGA.TE_TM(x0, 0, pol) << "\n";
				}

				term = template_funcs::DSQR(WGA.TE_TM(x0, 0, pol));

				sum += term;

				// update position
				x0 += dx;
			}

			if (loud) {
				write.close();
			}

			// compute last term in the sum
			last = template_funcs::DSQR(WGA.TE_TM(x0, 0, pol));

			integral = 0.5 * dx * (first + last + 2.0 * sum);

			if(scale) integral *= de_A; // multiply \Delta\epsilon_{a}

			integral *= 2.0;

			return integral;
		}
		else {
			return 0;

			std::string reason = "Error: double coupled_slabs::integrate_KAA(double pitch)\n";
			if (!wg_defined) reason += "WG parameters are not defined\n";
			//if (!c1) reason += "Integration Limits are not correct\n";
			if (!c3) reason += "WG pitch = " + template_funcs::toString(pitch, 3) + " is too small\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double coupled_slabs::integrate_KBB(double pitch, bool scale, bool loud)
{
	// compute the coupling integral of E_{by} E_{by}
	// range of integration is -a <= x <= +a

	try {
		bool c3 = pitch > 0.5 * (WA + WB) ? true : false;

		if (c3 && wg_defined) {
			int N = 501;
			double Lx = 2.5 * (WA + pitch + WB); // total length of simulation region
			double dx = Lx / (N - 1);
			double x0 = -0.5*WA; // start computing solutions here
			double xend = 0.5*WA; // stop computing solutions here
			double first, last, term, integral, sum = 0.0;

			std::string filename = "integrate_KBB_field_values.txt";
			std::ofstream write;

			if (loud) {
				write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);
			}

			// Evaluate the integral using the trapezoidal rule

			// compute first term in the sum
			first = template_funcs::DSQR(WGB.TE_TM(x0, 0, pol));

			// update position
			x0 += dx;

			// sum over the middle terms 
			while(x0 < xend) {

				if (loud) {
					write << std::setprecision(10) << x0 << " , " << WGB.TE_TM(x0-pitch, 0, pol) << " , " << WGB.TE_TM(x0 - pitch, 0, pol) << "\n";
				}

				term = template_funcs::DSQR(WGB.TE_TM(x0 - pitch, 0, pol));

				sum += term;

				// update position
				x0 += dx;
			}

			if (loud) {
				write.close();
			}

			// compute last term in the sum
			last = template_funcs::DSQR(WGB.TE_TM(x0, 0, pol));

			integral = 0.5 * dx * (first + last + 2.0 * sum);

			if(scale) integral *= de_B; // multiply \Delta\epsilon_{b}

			integral *= 2.0;

			return integral;
		}
		else {
			return 0;

			std::string reason = "Error: double coupled_slabs::integrate_KBB(double pitch)\n";
			if (!wg_defined) reason += "WG parameters are not defined\n";
			//if (!c1) reason += "Integration Limits are not correct\n";
			if (!c3) reason += "WG pitch = " + template_funcs::toString(pitch, 3) + " is too small\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double coupled_slabs::integrate_KAB(double pitch, bool scale, bool loud)
{
	// compute the coupling integral of E_{ay} E_{by}
	// range of integration is -a <= x <= +a

	try {
		bool c3 = pitch > 0.5 * (WA + WB) ? true : false;

		if (c3 && wg_defined) {
			int N = 501;
			double Lx = 2.5 * (WA + pitch + WB); // total length of simulation region
			double dx = Lx / (N - 1);
			double x0 = -0.5*WA; // start computing solutions here
			double xend = 0.5*WA; // stop computing solutions here
			double xoff = x0 - pitch; // use the offset value to evaluate mode solution in WGB
			double first, last, term, integral, sum = 0.0;

			std::string filename = "integrate_KAB_field_values.txt";
			std::ofstream write;

			if (loud) {
				write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);
			}

			// Evaluate the integral using the trapezoidal rule

			// compute first term in the sum
			first = (WGA.TE_TM(x0, 0, pol) * WGB.TE_TM(x0 - pitch, 0, pol));

			// update position
			x0 += dx;
			xoff += dx;

			// sum over the middle terms 
			while(x0 < xend) {

				if (loud) {
					write << std::setprecision(10) << x0 << " , " << WGA.TE_TM(x0, 0, pol) << " , " << WGB.TE_TM(x0-pitch, 0, pol) << "\n";
				}

				term = (WGA.TE_TM(x0, 0, pol) * WGB.TE_TM(x0 - pitch, 0, pol));

				sum += term;

				// update position
				x0 += dx;
				xoff += dx;
			}

			if (loud) {
				write.close();
			}

			// compute last term in the sum
			last = (WGA.TE_TM(x0, 0, pol) * WGB.TE_TM(x0 - pitch, 0, pol));

			integral = 0.5 * dx * (first + last + 2.0 * sum);

			if(scale) integral *= de_B; // multiply \Delta\epsilon_{b}

			integral *= 2.0;

			return integral;
		}
		else {
			return 0;

			std::string reason = "Error: double coupled_slabs::integrate_KAB(double pitch)\n";
			if (!wg_defined) reason += "WG parameters are not defined\n";
			//if (!c1) reason += "Integration Limits are not correct\n";
			if (!c3) reason += "WG pitch = " + template_funcs::toString(pitch, 3) + " is too small\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double coupled_slabs::integrate_KBA(double pitch, bool scale, bool loud)
{
	// compute the coupling integral of E_{ay} E_{by}
	// range of integration is D-b <= x <= D+b

	try {
		bool c3 = pitch > 0.5 * (WA + WB) ? true : false;

		if (c3 && wg_defined) {
			int N = 501;
			double Lx = 2.5 * (WA + pitch + WB); // total length of simulation region
			double dx = Lx / (N - 1);
			double x0 = pitch - 0.5 * WB; // start computing solutions here
			double xend = pitch + 0.5 * WB; // stop computing solutions here
			double xoff = x0 - pitch; // use the offset value to evaluate mode solution in WGB
			double first, last, term, integral, sum = 0.0;

			std::string filename = "integrate_KBA_field_values.txt";
			std::ofstream write;

			if (loud) {
				write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);
			}

			// Evaluate the integral using the trapezoidal rule

			// compute first term in the sum
			first = (WGA.TE_TM(x0, 0, pol) * WGB.TE_TM(x0 - pitch, 0, pol));

			// update position
			x0 += dx;
			xoff += dx;

			// sum over the middle terms 
			while (x0 < xend) {

				if (loud) {
					write << std::setprecision(10) << x0 << " , " << WGA.TE_TM(x0, 0, pol) << " , " << WGB.TE_TM(x0 - pitch, 0, pol) << "\n";
				}

				term = (WGA.TE_TM(x0, 0, pol) * WGB.TE_TM(x0 - pitch, 0, pol));

				sum += term;

				// update position
				x0 += dx;
				xoff += dx;
			}

			if (loud) {
				write.close();
			}

			// compute last term in the sum
			last = (WGA.TE_TM(x0, 0, pol) * WGB.TE_TM(x0 - pitch, 0, pol));

			integral = 0.5 * dx * (first + last + 2.0 * sum);

			if(scale) integral *= de_A; // multiply \Delta\epsilon_{a}

			integral *= 2.0;

			return integral;
		}
		else {
			return 0;

			std::string reason = "Error: double coupled_slabs::integrate_KBA(double pitch)\n";
			if (!wg_defined) reason += "WG parameters are not defined\n";
			//if (!c1) reason += "Integration Limits are not correct\n";
			if (!c3) reason += "WG pitch = " + template_funcs::toString(pitch, 3) + " is too small\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void coupled_slabs::define_P(double z, bool loud)
{
	// Compute the propagation matrix P at position z along the direction of propagation
	
	try {
	
		if (coeffs_defined) {

			P = vecut::zero_cmat(msize, msize);

			P[0][0] = exp(eye * bp * z);
			P[1][1] = exp(eye * bm * z);

			if (loud) {
				std::cout << "Propagation Matrix M z = " << z << "\n";
				vecut::print_cmat(P); 
				std::cout << "\n"; 
			}
		}
		else {
			std::string reason = "Error: void coupled_slabs::define_P(double z)\n";
			if (!coeffs_defined) reason += "propagation parameters are not defined\n";
			
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void coupled_slabs::define_M(double z, bool loud)
{
	// Compute the propagation matrix M = V . P . V^{-1} at position z along the direction of propagation

	try {

		if (coeffs_defined) {

			// define the propagation matrix P
			define_P(z, loud); 

			// Need to compute the matrix product here
			std::vector<std::vector<std::complex<double>>> M1; 
			
			M1 = vecut::cmat_cmat_product(V, P); 

			M = vecut::cmat_cmat_product(M1, Vinv); 

			if (loud) {
				std::cout << "Propagation Matrix M z = "<<z<<"\n";
				vecut::print_cmat(M);
				std::cout << "\n";
			}
		}
		else {
			std::string reason = "Error: void coupled_slabs::define_P(double z)\n";
			if (!coeffs_defined) reason += "propagation parameters are not defined\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void coupled_slabs::propagate(double length, double step_size, double a0, double b0, bool loud)
{
	// compute the field in the coupled waveguide along its length
	// assume initial condition (a0 b0)^{T}

	try {
		bool c1 = length > 0.0 ? true : false; 
		bool c4 = step_size > 0.0 && step_size <= length ? true : false; 
		bool c2 = a0 >= 0.0 && a0<= 1.0 ? true : false; 
		bool c3 = b0 >= 0.0 && b0 <= 1.0 ? true : false; 
		bool c10 = c1 && c2 && c3 && c4; 

		if (c10 && coeffs_defined) {

			std::vector<std::complex<double>> IC(msize, zero); 
			std::vector<std::complex<double>> AB(msize, zero); 

			IC[0] = a0; IC[1] = b0; 

			int Nsteps = 1 + static_cast<int>(length / step_size); 

			double z = 0; 

			std::string filename = "Coupled_Mode_Coefficients.txt";
			std::ofstream write;

			write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);

			for (int i = 0; i < Nsteps; i++) {
				
				define_M(z, loud); 

				// compute (a(z) b(z))^{T} = M.(a0 b0)^{T}
				AB = vecut::cmat_cvec_product(M, IC); 

				if (loud) {
					std::cout << "(a(z) b(z))^{T} z = " << z << "\n"; 
					for (int i = 0; i < msize; i++) std::cout << abs(AB[i]) << "\n"; 
					std::cout << "\n\n"; 
				}

				write << std::setprecision(10) << z << " , " << abs(AB[0]) << " , " << abs(AB[1]) << "\n"; 

				// compute the field profile in the waveguide
				if (i%10 == 0 && !Afield.empty() && !Bfield.empty() && Afield.size() == Bfield.size()) {

					std::string fieldfile = "Coupled_Field_dz_" + template_funcs::toString(step_size) + "_step_" + template_funcs::toString(i) + dottxt; 
					std::ofstream writefield; 

					writefield.open(fieldfile.c_str(), std::ios_base::out | std::ios_base::trunc);

					for (size_t j = 0; j < Afield.size(); j++) {
						writefield << std::setprecision(10) << pos[j] << " , " << abs( Afield[j]*AB[0] + Bfield[j]*AB[1] ) << "\n";
					}

					writefield.close(); 
				}

				// output the results

				z += step_size; 
			}

			write.close(); 

		}
		else {
			std::string reason = "Error: void coupled_slabs::propagate(double length, double a0, double b0)\n";
			if (!c1) reason += "length = " + template_funcs::toString(length, 2) + " is not positive\n"; 
			if (!c4) reason += "step_size = " + template_funcs::toString(step_size, 2) + " is too large\n";
			if (!c2) reason += "IC a0 = " + template_funcs::toString(a0, 2) + " is not in range\n"; 
			if (!c3) reason += "IC b0 = " + template_funcs::toString(b0, 2) + " is not in range\n"; 
			if (!coeffs_defined) reason += "propagation parameters are not defined\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}