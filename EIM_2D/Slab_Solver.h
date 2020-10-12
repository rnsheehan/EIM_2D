#ifndef SLAB_SOLVER_H
#define SLAB_SOLVER_H

// Slab Waveguide Solver
// R. Sheehan 8 - 7 - 2008

// Updated 18 - 7 - 2014
// Extra functionality previously associated with this class has been deprecated
// This class will now only perform the calculation of slab waveguide propagation constants and mode profiles
// Other types of analysis will inherit this class in future
// R. Sheehan 18 - 7 - 2014

class slab_wg{
public:
	//Constructor
	slab_wg(void);
	slab_wg(double width, double lambda, double ncore, double nsub, double nclad);

	// Destructor	
	~slab_wg(); // This is called implicitly once the class has gone out of scope

	//Setters
	void set_params(double width, double lambda, double ncore, double nsub, double nclad);

	void clearbeta(bool t);//Remove all elements from the vector
	
	//Getters
	int nbeta(bool mode);

	// Getters
	double _beta(int i,bool t); // returns a particular propagation constant

	double _neff(int i,bool t); // returns the effective index

	//Member Functions
	void calculate_all_modes(int N, double Lx, std::string &storage_directory); // compute modes for both polarisations and output the results

	void calculate_mode(bool mode, int mode_num, int N, double Lx); // compute the modes for a single polarisation

	void calculate_neff(bool mode); // compute effective indices for one polarisation

	void calculate_all_neffs(std::string &storage_directory); // compute effective indices for both polarisations, output results

	void neff_search(bool mode); // Self-explanatory

	void output_modes(bool mode, int N, double Lx, std::string &storage_directory); // Output solutions to a file
	
	void output_all_stats(std::string &storage_directory); // Output calculation data to a file

private:
	// methods that the user does not need access to

	// Functions needed to compute the shape of the waveguide mode
	double h(int i,bool t); // wavenumber in core

	double p(int i,bool t); // wavenumber in substrate

	double q(int i,bool t); // wavenumber in cladding

	double g1(int i,bool t);

	double g2(int i,bool t);

	double deff(int i,bool t); // effective width of waveguide mode

	double phase(int i,bool t); // phase of waveguide mode

	double norm_const(int i,bool t); // normalisation constant of waveguide mode

	double conf_fact(int i,bool t); // Mode confinement factor
	
	double TE_TM(double x,int i,bool mode); // shape of waveguide mode

	double eigeneqn(double x,bool t); // Non-linear equation for the propagation constants

	double eigeneqn_3(double x,bool t,int mm); // Non-linear equation for the propagation constants
	
	double zbrent(double x1,double x2,double tol,bool t,int mm); // Brent method search for roots of eigeneqn_3

	void output_stats(bool mode, std::ofstream &file_obj); // write computed mode statistics to a file

private:
	int M; // Predicted number of modes that waveguide can support

	double d; // Waveguide Width
	double nc; // Core Index
	double ncl; // Cladding Index
	double ns; // Substrate Index
	double g; // Asymmetry factor	

	double aa; // ratio of (nc / ns)^{2}
	double bb; // ratio of (nc / ncl)^{2}

	double k_sqr; // k_{0}^{2}
	double nc_sqr; // n_{c}^{2}
	double ns_sqr; // n_{s}^{2}
	double ncl_sqr; // n_{cl}^{2}

	double k_sqr_nc_sqr; // k_{0}^{2} n_{c}^{2}
	double k_sqr_ns_sqr; // k_{0}^{2} n_{s}^{2}
	double k_sqr_ncl_sqr; // k_{0}^{2} n_{cl}^{2}
	
	double l; // Wavelength
	double V; // V-parameter
	double na; // Numerical Aperture
	double k; // wavenumber
	double w; // frequency
	double lower; // Lower bound on the propagation constant
	double upper; // Upper bound on the propagation constant
	double efieldint; // Coefficient of the intensity of the electric field
	double hfieldint; // Coefficienct of the intensity of the magnetic field
		
	std::vector<double> betaE; // TE propagation constants
	std::vector<double> betaH; // TM propagation constants
};

#endif