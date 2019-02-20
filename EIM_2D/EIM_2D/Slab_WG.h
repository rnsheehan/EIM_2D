#ifndef SLAB_WG_H
#define SLAB_WG_H

// Declaration of the classes needed to make the slab waveguide solver
// R. Sheehan 2 - 5 - 2018

class slab {
public:
	// Constructors
	slab(void);

	// Deconstructor
	~slab(void);

protected:	
	int nbeta(bool mode); // get the number of computed modes in a waveguide

	double h(int i, bool t); // wavenumber in core

	double p(int i, bool t); // wavenumber in substrate

	double q(int i, bool t); // wavenumber in cladding

	double r(int i, bool t); // wavenumber in rib

	double _beta(int i, bool t); // return the computed propagation constants

protected:
	// protected members, only to be accessed through derived classes
	int M; // Predicted number of modes that waveguide can support

	double d; // Waveguide Width
	double dr; // Rib Width
	double nc; // Core Index
	double ncl; // Cladding Index
	double ns; // Substrate Index
	double nr; // Rib Index, only used by Four layer slab
	double nm; // nm = Max(nr,ns)

	double l; // Wavelength
	double V; // V-parameter
	double na; // Numerical Aperture
	double k; // wavenumber
	double w; // frequency

	double lower; // Lower bound on the propagation constant
	double upper; // Upper bound on the propagation constant

	double efieldint; // Coefficient of the intensity of the electric field
	double hfieldint; // Coefficienct of the intensity of the magnetic field

	double k_sqr; // k_{0}^{2}
	double nc_sqr; // n_{c}^{2}
	double ns_sqr; // n_{s}^{2}
	double ncl_sqr; // n_{cl}^{2}
	double nr_sqr; // n_{r}^{2}
	double nm_sqr; // n_{m}^{2}

	double aa; // ratio of (nc / ns)^{2}
	double bb; // ratio of (nc / ncl)^{2}

	double etacs; // n_{c}^{2} / n_{s}^{2}
	double etacr; // n_{c}^{2} / n_{r}^{2}
	double etarcl; // n_{r}^{2} / n_{cl}^{2}

	double k_sqr_nc_sqr; // k_{0}^{2} n_{c}^{2}
	double k_sqr_ns_sqr; // k_{0}^{2} n_{s}^{2}
	double k_sqr_ncl_sqr; // k_{0}^{2} n_{cl}^{2}
	double k_sqr_nr_sqr; // k_{0}^{2} n_{r}^{2}
	double k_sqr_nm_sqr; // k_{0}^{2} n_{m}^{2}

	std::vector<double> betaE; // TE propagation constants
	std::vector<double> betaH; // TM propagation constants
};

// Three Layer Slab

class slab_tl_neff : protected slab {
	// class which is used to compute the effective indices of a three layer slab waveguide
public:
	// Constructor
	slab_tl_neff(void); 

	slab_tl_neff(double width, double lambda, double ncore, double nsub, double nclad);

	void set_params(double width, double lambda, double ncore, double nsub, double nclad); // assign the parameters needed to perform the neff_search calculation

	void neff_search(bool mode); // compute the effective indices for a waveguide

private:
	// methods that the user does not need access to

	double eigeneqn_3(double x, bool t, int mm); // Non-linear equation for the propagation constants

	double zbrent(double x1, double x2, double tol, bool t, int mm); // Brent method search for roots of eigeneqn_3

private:
	double g; // Asymmetry factor
};

class slab_tl_mode : public slab_tl_neff {
	// class which is used to compute the shapes of optical modes in a three layer slab
public:
	slab_tl_mode(void); 

	slab_tl_mode(double width, double lambda, double ncore, double nsub, double nclad);

	void compute_neff(bool mode);

	void output_all_stats(std::string &storage_directory); 

	void output_modes(bool mode, int N, double Lx, std::string &storage_directory); // Output solutions to a file

private:
	// Functions needed to compute the shape of the waveguide mode
	double g1(int i, bool t);

	double g2(int i, bool t);

	double deff(int i, bool t); // effective width of waveguide mode

	double phase(int i, bool t); // phase of waveguide mode

	double norm_const(int i, bool t); // normalisation constant of waveguide mode

	double conf_fact(int i, bool t); // Mode confinement factor

	double TE_TM(double x, int i, bool mode); // shape of waveguide mode

	double eigeneqn(double x, bool t); // Non-linear equation for the propagation constants
	
	void output_stats(bool mode, std::ofstream &file_obj); // write computed mode statistics to a file
};

// Four Layer Slab

class slab_fl_neff_B : protected slab {
	// class which is used to compute the effective indices in a Case B four layer slab

	// Case A => Field Oscillating in Core and Ridge
	// For there to be a solution one has to have ns <= ncl < nr < nc
	// Case A is not included in this project
	
	// Case B: Field Oscillating in Core Only
	// For there to be a solution one has to have ncl < nm < nc, where nm = Max(nr,ns)
public:
	// Constructors
	slab_fl_neff_B(void);

	slab_fl_neff_B(double width, double rib_width, double lambda, double ncore, double nsub, double nclad, double nrib);

	// Methods
	void set_params(double width, double rib_width, double lambda, double ncore, double nsub, double nclad, double nrib);

	void neff_search(bool mode);

private:
	double eigeneqn_b(double x, int mm, bool t);

	double zbrent(double x1, double x2, double tol, bool t, int mm); // Brent method search for roots of eigeneqn_b
};

class slab_fl_mode_B : public slab_fl_neff_B {
	// Class for computing the effective indices and mode profiles in type B four layer slab
public:
	slab_fl_mode_B(void);

	slab_fl_mode_B(double width, double rib_width, double lambda, double ncore, double nsub, double nclad, double nrib);

	void compute_neff(bool mode); 

	//void output_all_stats(std::string &storage_directory);

	//void output_modes(bool mode, int N, double Lx, std::string &storage_directory); // Output solutions to a file

private:
	double phase(int i, bool t); // phase of waveguide mode

	double TE_TM(double x, int i, bool mode); // shape of waveguide mode

	//double eigeneqn(double x, bool t); // Non-linear equation for the propagation constants

	//void output_stats(bool mode, std::ofstream &file_obj); // write computed mode statistics to a file
};

// class for computing the coupling coefficient of the coupled slab waveguide

class coupled_slab_tl_neff : public slab_tl_neff {

public:
	coupled_slab_tl_neff();

	coupled_slab_tl_neff(double separation, double width, double lambda, double ncore, double nsub);

	void set_params(double separation, double width, double lambda, double ncore, double nsub);

	double compute_coupling_coeff(bool mode); 

private:
	double slab_sep; 
	double coupling_coeff; 
	double L_coupling; 
};

#endif