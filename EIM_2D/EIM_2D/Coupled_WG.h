#ifndef COUPLED_WG_H
#define COUPLED_WG_H

// Declaration of class to be used for performing approximate analysis of coupled 2D waveguides
// The calculations reduces the 2D waveguides to 1D slab using the EIM
// Coupling calculations are then carried out in 1D
// The idea is to approximate the coupling coefficient between the waveguides
// R. Sheehan 4 - 12 - 2020

// Dev. Note
// A constraint of the coupled-slab calculation is that the slabs must share the same substrate and that they must be symmetric
// Now that I think about it this is overly constraining the coupled slabs need not be symmetric
// The coupled slabs must have common material between them, but their outer materials can be different
// i.e. it should be possible to do the coupled slab calculations with non-symmetric slabs as long as the 
// cladding RI of WGA matches the substrate RI of WGB, a picture will make more sense here but I'm fairly sure I'm right. 
// if you extend the coupled slab calculation to non-symmetric waveguides then you have the possibility of looking at VERTICAL coupling as well as lateral coupling
// vertical coupling between waveguides is a very interesting topic as you know. 
// The advantage of studying vertical coupling using this approach is that it is quasi-analytical and very-fast
// Worthty of further investigation given the applications
// R. Sheehan 4 - 12 - 2020

// Dev. Note
// It might be preferable to set up derived classes for the different cases of 
// 1. waveguide coupling to a copy of itself
// 2. waveguide coupling to a different waveguide
// In case 1 you could reduce the overhead by reducing the number of slab waveguide mode calculations
// Is it worth the effort though? The overhead due to the calculation is tiny
// Even then is the effort needed to define the derived classes worth it
// Just trying to remember the logic around virtual functions is making my head hurt
// You can always just run the calculation by passing two pointers to the same waveguide object
// This is in fact what I have done in the "derived class"
// R. Sheehan 4 - 12 - 2020

class coupWG {
public:
	coupWG();
	coupWG(EIM* wgobjA, EIM* wgobjB); 
	coupWG(coupWG& obj); 
	~coupWG();

	void set_params(EIM* wgobjA, EIM* wgobjB);

	void reduce_wg(bool loud = false); 

	void coupling_coeffs(double pitch, bool loud = false);

protected:
	bool params_defined; 
	bool waveguides_reduced; 

	EIM* WGA; // pointer to 2D waveguide A
	EIM* WGB; // pointer to 2D waveguide B

	coupled_slabs coupling; // object used to perform the coupling calculations
};

// Derived class that can be used to examine the coupling characterisics of a waveguide with a copy of itself
// The only difference is that the constructor passes the same pointer wgobj to both arguments of set_params
// I guess it makes a slightly less complicated interface available, but it's just a wrapper
// R. Sheehan 4 - 12 - 2020

class coupWGitself : public coupWG {
public:
	coupWGitself();
	coupWGitself(EIM* wgobj);
};

#endif 

