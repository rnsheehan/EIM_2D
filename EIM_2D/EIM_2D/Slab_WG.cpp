#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definitions for the slab base class

slab::slab(void)
{
	// Default Constructor
	M = 0;

	k_sqr = nc_sqr = ncl_sqr = ns_sqr = nr_sqr = nm_sqr = k_sqr_nc_sqr = k_sqr_ns_sqr = k_sqr_ncl_sqr = k_sqr_nr_sqr = 0.0;
	k_sqr_nm_sqr = aa = bb = d = dr = l = nc = ns = ncl = nr = nm = etacr = etacs = etarcl = 0.0;
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

double slab::r(int i, bool t)
{
	//Wavenumber in Rib

	try {

		if (nbeta(t) > 0) {

			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab::r(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt( template_funcs::DSQR(betaE[i]) - k_sqr_nr_sqr);
				}
				else {
					return sqrt( template_funcs::DSQR(betaH[i]) - k_sqr_nr_sqr);
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

	for (int i = 0; i<nbeta(mode); i++) {
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
				throw std::range_error("Error: double slab_wg::get_neff(int i, bool t)\n Attempting to access arrays out of range\n");
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
			throw std::invalid_argument("Error: double slab::get_neff(int i, bool t)\nNo modes have been computed\n");
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
		double db = ((upper - lower) / (100 - 1));

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
				double t = d / 2;

				if (x > t) {
					return norm_const(i, mode)*cos(h(i, mode)*t - phase(i, mode))*exp(-q(i, mode)*(x - t)); // Cladding
				}
				else if (x < -t) {
					return norm_const(i, mode)*cos(h(i, mode)*t + phase(i, mode))*exp(p(i, mode)*(x + t)); // Substrate
				}
				else {
					return norm_const(i, mode)*cos(h(i, mode)*x - phase(i, mode)); // Core
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl_mode::TE_TM(int i, bool t)\nNo modes have been computed\n");
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
			M = static_cast<int>(std::max(1.0, ceil((2.0*V / PI))));

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
			//if (!t) q *= etarcl; // multiply q by etarcl in the case of TM polarisation

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
			double fa = eigeneqn_b(a, t, mm), fb = eigeneqn_b(b, t, mm);

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
				fb = eigeneqn_b(b, t, mm);
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
					double rr = r(i, t); 
					return -1.0*( ( dr * rr ) + atanh(qq/rr) );
				}
				else {
					double qq = q(i, t);
					double rr = r(i, t);
					return -1.0*( ( dr * rr ) + atanh( (etarcl * qq) / rr) );
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
				double rr = r(i, mode); 
				double wh = d * hh;
				double qq = q(i, mode);
				double pp = p(i, mode); 
				double ph = phase(i, mode); 
				double rrhh = mode ? (rr / hh) * tan(ph) : (rr / hh) * tan(ph) * etacr;
				if (x < -d) {
					// x < -W
					return ( ( cos(wh) + ( rrhh * sin(wh) ) )*exp( pp * ( x + d ) ) );
				}
				else if (x >= -d && x <= 0) {
					// -W < x < 0
					return ( cos( hh * x ) - rrhh * sin( hh * x ) ); 
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
	// returned value is in units of um^{-1}

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