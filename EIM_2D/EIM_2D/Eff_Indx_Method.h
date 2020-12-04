#ifndef EFF_INDX_METHOD_H
#define EFF_INDX_METHOD_H

// Implementation of the effective index for 2D waveguides
// Uses three and four layer slab calculations
// This implementation is a refactoring of an older scheme that has become too cumbersome and unwieldy
// R. Sheehan 20 - 2 - 2019

// EIM is the base class
// EIM will contain all parameters for the derived classes as well as the implementation of the get_index which will be required
// by all derived classes. Each derived class will have an instance of the function reduce_wg which will reduce the 2D WG structure
// to an equivalent three layer slab structure

// You should be able to reduce the number of set_params methods for the derived classes to a single method defined in EIM base class
// Although it may be better to have them defined explicitly for the derived classes so that each calculation in reduce_wg is correctly
// defined. Otherwise you need to have special cases for each of the different waveguides, and that defeats the purpose of using dynamic binding
// Sticking with the multiple set_params methods!

// What's the best way to implement the EIM for coupled waveguides? 
// Should you define EIM_Coupled as a derived class of EIM and then subsequently define 
// Coupled_Rect as derived class of Rect? Or should you just go for Coupled_Rect as 
// derived class of Rect with members Rect WG1 and Rect WG2
// Coupling calculations can be done by having an instance of Coupled_Slabs
// Might be best just to work out a few examples manually and see if a more generic structure
// can be deduced from that. This is probably the easiest way to do it tbh. 
// R. Sheehan 2 - 12 - 2020

// There's a better way. See testing::eim_arb_wg()
// You can set up an object that takes two pointers to EIM objects as parameters
// The EIM object is a pointer to a derived  waveguide object (Rect, Rib, Ridge etc)
// The calculations needed to reduce the 2D WG to a 1D slab can be performed by the EIM object pointer
// This will work for any pair of pointer to waveguide objects so no need to brute force it for each case
// It's then very straightforward to use coupled_slab_tl_neff to perform the coupled waveguide analysis
// All thanks to object oriented programming
// R. Sheehan 3 - 12 - 2020

class wg_dims {
public:
	wg_dims();
	wg_dims(wg_dims &obj); 

	void set_dims(wg_dims &obj); 

	void set_rect_wire(double W, double H);
	void set_rib(double W, double E, double T);
	void set_ridge(double W, double E, double T, double D);

	// getters
	inline bool defined() { return params_defined; }

	inline int get_wg_code() { return wg_type;  }

	inline double get_W() { return width;  }
	inline double get_H() { return height; }
	inline double get_E() { return etch_depth; }
	inline double get_T() { return slab_height; }
	inline double get_D() { return core_height; }

private:
	bool params_defined; 

	int wg_type; // code to identify the waveguide geometry

	// Waveguide Dimensions
	double width; // W rib / ridge width in units of um
	double height; // H total height in units of um
	double etch_depth; // E rib / ridge etch depth in units of um
	double slab_height; // T rib / ridge slab height in units of um
	double core_height; // D ridge core layer thickness in units of um
};

class ri_vals {
public: 
	ri_vals();

	ri_vals(ri_vals &obj); 

	void set_ri(ri_vals &obj);

	void set_rect(double Ncore, double Nclad, double WL);
	void set_rib_wire(double Ncore, double Nsub, double Nclad, double WL);
	void set_ridge(double Ncore, double Nsub, double Nrib, double Nclad, double WL);

	// getters
	inline bool defined() { return params_defined;  }
	inline double get_Ncore() { return ncore; }
	inline double get_Nsub() { return nsub; }
	inline double get_Nrib() { return nrib; }
	inline double get_Nclad() { return nclad; }
	inline double get_WL() { return lambda;  }

private:
	bool params_defined; 
	// Waveguide refractive indices
	double ncore; // core RI
	double nsub; // substrate RI
	double nrib; // rib RI
	double nclad; // claddin RI

	double lambda; // Wavelength in units of um
};

class EIM
{
public:
	EIM();
	~EIM(); 

	// By right you should have another virtual function set_params, however, this would be awkward as hell to implement 
	// since each derived class requires a different number of parameters, best just to have set_params declared and defined
	// in each derived class

	virtual void reduce_wg() = 0; // each derived class will reduce 2D - 1D slightly differently
	virtual void set_params() = 0; 

	virtual void set_params(bool polarisation, wg_dims &, ri_vals &);

	void get_index(bool loud); // compute the effective index of the 2D structure

	double neff_value(); // return the effective index of the fundamental mode

	inline bool defined() { return params_defined;  }
	inline std::vector<double> reduced_cladding() { return eta_one; }
	inline std::vector<double> reduced_core() { return eta_two; }
	inline wg_dims get_dims() { return dimensions; }
	inline ri_vals get_RI() { return ref_indx; }

	//double coupling_coefficient(double &separ); // estimate the coupling coefficient between this waveguide and a copy of itself, possibly best to deprecate this

	// setters
	/*inline void set_W(double &val) { width = val;  }
	inline void set_H(double &val) { height = val; }
	inline void set_E(double &val) { etch_depth = val; }
	inline void set_T(double &val) { slab_height = val; }
	inline void set_D(double &val) { core_height = val; }

	inline void set_Ncore(double &val) { ncore = val; }
	inline void set_Nsub(double &val) { nsub = val; }
	inline void set_Nrib(double &val) { nrib = val; }
	inline void set_Nclad(double &val) { nclad = val; }
	
	inline void set_WL(double &val) { lambda = val; }*/

protected:
	bool pol; // polarisation
	bool not_pol; // !polarisation
	bool params_defined; 

	// Waveguide Dimensions
	double width; // W rib / ridge width in units of um
	double height; // H total height in units of um
	double etch_depth; // E rib / ridge etch depth in units of um
	double slab_height; // T rib / ridge slab height in units of um
	double core_height; // D ridge core layer thickness in units of um

	// Waveguide refractive indices
	double ncore; // core RI
	double nsub; // substrate RI
	double nrib; // rib RI
	double nclad; // cladding RI

	double lambda; // Wavelength in units of um

	wg_dims dimensions; // Object for accessing all waveguide dimension parameters

	ri_vals ref_indx; // Object for accessing all waveguide refractive index parameters

	// vectors to hold the reduced index values for the side-slab and the core region
	// reduced waveguide is symmetric slab whose cladding is RI eta_one and core RI is eta_two
	std::vector<double> eta_one; // reduced cladding indices
	std::vector<double> eta_two; // reduced core indices

	slab_tl_neff TLS; // three layer slab object

	slab_fl_neff_B FLS; // four layer slab object
};

class Rectangular : public EIM 
{
public:
	Rectangular(); 

	Rectangular(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj);

	void set_params(); 
	//void set_params(bool polarisation, double W, double H, double Ncore, double Nclad, double WL);
	void set_params(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj);
 
	void reduce_wg(); 
};

class Wire : public EIM
{
public:
	Wire();

	Wire(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj);

	void set_params();
	//void set_params(bool polarisation, double W, double H, double Ncore, double Nsub, double Nclad, double WL);
	void set_params(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj);

	void reduce_wg();
};

class Rib : public EIM
{
public:
	Rib(); 

	Rib(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj);

	void set_params();
	//void set_params(bool polarisation, double W, double E, double T, double Ncore, double Nsub, double Nclad, double WL);
	void set_params(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj);

	void reduce_wg(); 
};

class Shallow_Ridge : public EIM
{
public:
	Shallow_Ridge();

	Shallow_Ridge(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj);

	void set_params();
	//void set_params(bool polarisation, double W, double E, double T, double D, double Ncore, double Nsub, double Nrib, double Nclad, double WL);
	void set_params(bool polarisation, wg_dims &dim_obj, ri_vals &ri_obj);

	void reduce_wg();
};

#endif
