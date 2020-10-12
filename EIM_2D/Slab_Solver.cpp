#ifndef ATTACH_H
#include "Attach.h"
#endif

// Slab Waveguide Solver
// R. Sheehan 8 - 7 - 2008

// Updated R. Sheehan 31 - 8 - 2016

//Constructor
slab_wg::slab_wg()
{
	// Default constructor

	M = 0; 

	k_sqr = nc_sqr = ncl_sqr = ns_sqr = k_sqr_nc_sqr = k_sqr_ns_sqr = k_sqr_ncl_sqr = 0.0; 
	aa = bb = d = l = nc = ns = ncl = 0.0; 
	g = V = na = k = lower = upper = upper = w = efieldint = hfieldint = 0.0; 
}

slab_wg::slab_wg(double width, double lambda, double ncore, double nsub, double nclad)
{
	// Primary constructor

	set_params(width, lambda, ncore, nsub, nclad); 
}

slab_wg::~slab_wg()
{
	// Deconstructor

	betaE.clear(); 

	betaH.clear(); 
}

void slab_wg::set_params(double width,double lambda,double ncore,double nsub,double nclad)
{ 
	try{

		bool c1 = width > 0.0 ? true : false; 
		bool c2 = lambda > 0.0 ? true : false; 
		bool c3 = nclad >= 1.0 ? true : false; 
		bool c4 = nsub >= 1.0 ? true : false; 
		bool c5 = ncore > std::max(nsub, nclad) ? true : false; 

		if(c1 && c2 && c3 && c4 && c5){

			d=width;

			l = lambda;

			nc = ncore;
			nc_sqr = template_funcs::DSQR(nc); 

			ns = std::max(nsub, nclad);
			ns_sqr = template_funcs::DSQR(ns); 
			
			ncl = std::min(nsub, nclad);
			ncl_sqr = template_funcs::DSQR(ncl); 

			aa = ( nc_sqr / ns_sqr ); 

			bb = ( nc_sqr / ncl_sqr );

			k = Two_PI/l; 
			k_sqr = template_funcs::DSQR(k); // k_{0}^{2}

			k_sqr_nc_sqr = k_sqr*nc_sqr; // k_{0}^{2} n_{c}^{2}
			k_sqr_ns_sqr = k_sqr*ns_sqr; // k_{0}^{2} n_{c}^{2}
			k_sqr_ncl_sqr = k_sqr*ncl_sqr; // k_{0}^{2} n_{c}^{2}

			double x = nc_sqr - ns_sqr; 
			double y = ns_sqr - ncl_sqr; 

			// WG asymmetry paramater
			// g = ( template_funcs::DSQR(ns) - template_funcs::DSQR(ncl) )/( template_funcs::DSQR(nc) - template_funcs::DSQR(ns) );
			g = y > 0.0 ? ( y / x ) : 0.0 ; 

			// na = sqrt( template_funcs::DSQR(nc) - template_funcs::DSQR(ns) );
			na = sqrt(x); // numerical aperture
	
			//V = ( PI*d*sqrt( template_funcs::DSQR(nc) - template_funcs::DSQR(ns) ) ) / l;
			V = ( PI*d*na )/l; // V-parameter
 
			// predicted number of modes
			M = y > 0 ? static_cast<int>( std::max( 1.0, ceil( (2.0*V/PI) - atan(g)/PI ) ) ) : static_cast<int>( std::max( 1.0, ceil( (2.0*V/PI) ) ) );

			lower = k*ns; // lower bound of search space

			upper = k*nc; // upper bound of search space

			w = k*SPEED_OF_LIGHT;

			efieldint = 0.5*EPSILON*SPEED_OF_LIGHT;

			hfieldint = 0.5*MU*SPEED_OF_LIGHT;

			//Empty the std::vector each time a new instance of the class is called
			betaE.clear();
			betaH.clear();		
		}
		else{
			std::string reason = "Error: void slab_wg::set_params(double width,double lambda,double ncore,double nsub,double nclad)\n"; 
			if(!c1) reason += "WG width = " + template_funcs::toString(width, 3) + " is negative\n"; 
			if(!c2) reason += "Wavelength = " + template_funcs::toString(lambda, 3) + " is negative\n";
			if(!c3) reason += "Cladding Index = " + template_funcs::toString(nclad, 3) + " is less than one\n";
			if(!c4) reason += "Substrate Index = " + template_funcs::toString(nsub, 3) + " is less than one\n";
			if(!c5) reason += "Core Index = " + template_funcs::toString(ncore, 3) + " is less than cladding / substrate index\n";

			throw std::invalid_argument(reason); 
		}	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void slab_wg::clearbeta(bool t)
{
	if(t){
		betaE.clear();
	}
	else{
		betaH.clear();
	}
}

//Getters
int slab_wg::nbeta(bool mode)
{
	// return the number of computed propagation constants for a given polarisation
	// it can be the case that the predicted value M is greater than the actual number of modes in a waveguide

	return static_cast<int>( mode ? betaE.size() : betaH.size() ); 

	/*if(mode){
		return static_cast<int>(betaE.size());
	}
	else{
		return static_cast<int>(betaH.size());
	}*/
}

double slab_wg::_beta(int i, bool t)
{
	// return the i^{th} propagation constant for a given polarisation
	// by convention t == true => TE modes, t == false => TM modes

	try{
		if(i<0 || i > nbeta(t) ){
			throw std::range_error("Error: double slab_wg::_beta(int i, bool t)\n Attempting to access arrays out of range\n"); 
			return 0;
		}
		else{
			if(t){
				return betaE[i];
			}
			else{
				return betaH[i];
			}
		}	
	}
	catch(std::range_error &e){
		std::cerr<<e.what();
		return 0; 
	}
}

double slab_wg::_neff(int i, bool t)
{
	// return the i^{th} effective index for a given polarisation
	// by convention t == true => TE modes, t == false => TM modes

	try{

		if(i<0 || i > nbeta(t)){
			throw std::range_error("Error: double slab_wg::_neff(int i, bool t)\n Attempting to access arrays out of range\n"); 
			return 0;
		}
		else{
			if(t){
				return betaE[i] / k;
			}
			else{
				return betaH[i] / k;
			}
		}	
	}
	catch(std::range_error &e){
		std::cerr<<e.what();
		return 0; 
	}
}

double slab_wg::h(int i,bool t)
{
	//Wavenumber in Core

	try{

		if(i<0 || i > nbeta(t)){
			throw std::range_error("Error: double slab_wg::h(int i, bool t)\n Attempting to access arrays out of range\n"); 
			return 0;
		}
		else{
			if(t){
				return sqrt( k_sqr_nc_sqr - template_funcs::DSQR(betaE[i]) );
			}
			else{
				return sqrt( k_sqr_nc_sqr - template_funcs::DSQR(betaH[i]) );
			}
		}	
	}
	catch(std::range_error &e){
		std::cerr<<e.what();
		return 0; 
	}
}

double slab_wg::p(int i,bool t)
{
	//Wavenumber in Substrate

	try{

		if(i<0 || i > nbeta(t)){
			throw std::range_error("Error: double slab_wg::p(int i, bool t)\n Attempting to access arrays out of range\n"); 
			return 0;
		}
		else{
			if(t){
				return sqrt( template_funcs::DSQR(betaE[i]) - k_sqr_ns_sqr);
			}
			else{
				return sqrt( template_funcs::DSQR(betaH[i]) - k_sqr_ns_sqr );
			}
		}	
	}
	catch(std::range_error &e){
		std::cerr<<e.what();
		return 0; 
	}
}

double slab_wg::q(int i,bool t)
{
	//Wavenumber in Cladding

	try{

		if(i<0 || i > nbeta(t)){
			throw std::range_error("Error: double slab_wg::q(int i, bool t)\n Attempting to access arrays out of range\n"); 
			return 0;
		}
		else{
			if(t){
				return sqrt( template_funcs::DSQR(betaE[i]) - k_sqr_ncl_sqr);
			}
			else{
				return sqrt( template_funcs::DSQR(betaH[i]) - k_sqr_ncl_sqr );
			}
		}	
	}
	catch(std::range_error &e){
		std::cerr<<e.what();
		return 0; 
	}
}

double slab_wg::g1(int i,bool t)
{
	try{

		if(i<0 || i > nbeta(t)){
			throw std::range_error("Error: double slab_wg::g1(int i, bool t)\n Attempting to access arrays out of range\n"); 
			return 0;
		}
		else{
			if(t){
				//return (template_funcs::DSQR(betaE[i])/(template_funcs::DSQR(k)*template_funcs::DSQR(nc)))+(template_funcs::DSQR(betaE[i])/(template_funcs::DSQR(k)*template_funcs::DSQR(ns)))-1;
				double beta_sqr = template_funcs::DSQR(betaE[i]); 
				return ( beta_sqr / k_sqr_nc_sqr ) + ( beta_sqr / k_sqr_ns_sqr ) - 1.0; 
			}
			else{
				//return (template_funcs::DSQR(betaH[i])/(template_funcs::DSQR(k)*template_funcs::DSQR(nc)))+(template_funcs::DSQR(betaH[i])/(template_funcs::DSQR(k)*template_funcs::DSQR(ns)))-1;
				double beta_sqr = template_funcs::DSQR(betaH[i]); 
				return ( beta_sqr / k_sqr_nc_sqr ) + ( beta_sqr / k_sqr_ns_sqr ) - 1.0; 
			}
		}	
	}
	catch(std::range_error &e){
		std::cerr<<e.what();
		return 0; 
	}
}

double slab_wg::g2(int i,bool t)
{
	try{
		
		if(i<0 || i > nbeta(t) ){
			throw std::range_error("Error: double slab_wg::g1(int i, bool t)\n Attempting to access arrays out of range\n"); 
			return 0;
		}
		else{
			if(t){
				//return (template_funcs::DSQR(betaE[i])/(template_funcs::DSQR(k)*template_funcs::DSQR(nc)))+(template_funcs::DSQR(betaE[i])/(template_funcs::DSQR(k)*template_funcs::DSQR(ncl)))-1;
				double beta_sqr = template_funcs::DSQR(betaE[i]); 
				return (  beta_sqr / k_sqr_nc_sqr ) + ( beta_sqr / k_sqr_ncl_sqr ) - 1.0; 
			}
			else{
				//return (template_funcs::DSQR(betaH[i])/(template_funcs::DSQR(k)*template_funcs::DSQR(nc)))+(template_funcs::DSQR(betaH[i])/(template_funcs::DSQR(k)*template_funcs::DSQR(ncl)))-1;
				double beta_sqr = template_funcs::DSQR(betaH[i]); 
				return ( beta_sqr / k_sqr_nc_sqr ) + ( beta_sqr / k_sqr_ncl_sqr ) - 1.0; 
			}
		}	
	}
	catch(std::range_error &e){
		std::cerr<<e.what();
		return 0; 
	}
}

double slab_wg::deff(int i,bool t)
{
	// Effective width of waveguide mode

	try{

		if(i<0 || i > nbeta(t)){
			throw std::range_error("Error: double slab_wg::deff(int i, bool t)\n Attempting to access arrays out of range\n"); 
			return 0;
		}
		else{
			if(t){
				return d + ( 1.0 / p(i,t) ) + ( 1.0 / q(i,t) );
			}
			else{
				return d + ( 1.0 / ( g1(i,t) * p(i,t) ) ) + ( 1.0 / ( g2(i,t) * q(i,t) ) );
			}
		}	
	}
	catch(std::range_error &e){
		std::cerr<<e.what();
		return 0; 
	}
}

double slab_wg::phase(int i,bool t)
{
	// Phase of mode in slab waveguide

	try{
		
		if(i<0 || i > nbeta(t)){
			throw std::range_error("Error: double slab_wg::phase(int i, bool t)\n Attempting to access arrays out of range\n"); 
			return 0;
		}
		else{
			if(t){
				double hh = h(i, t); 
				double pp = p(i, t); 
				double qq = q(i, t); 
				return i * (PI/2) + 0.5 * atan( ( pp / hh ) ) - 0.5 * atan( ( qq / hh ) );
			}
			else{
				//return i*(PI/2)+0.5*atan((template_funcs::DSQR(nc)/template_funcs::DSQR(ns))*(p(i,t)/h(i,t)))-0.5*atan((template_funcs::DSQR(nc)/template_funcs::DSQR(ncl))*(q(i,t)/h(i,t)));
				double hh = h(i, t); 
				double pp = p(i, t); 
				double qq = q(i, t); 
				return  i * (PI/2) + 0.5*atan( aa*( pp / hh ) ) - 0.5*atan( bb * ( qq / hh ) );
			}
		}	
	}
	catch(std::range_error &e){
		std::cerr<<e.what();
		return 0; 
	}
}

double slab_wg::norm_const(int i,bool t)
{
	// Normalisation constant for mode in a slab waveguide

	try{
		
		if(i<0 || i > nbeta(t)){
			throw std::range_error("Error: double slab_wg::norm_const(int i, bool t)\n Attempting to access arrays out of range\n"); 
			return 0;
		}
		else{
			if(t){
				return sqrt( ( w * MU ) / ( betaE[i] * deff(i,t) ) );
			}
			else{
				//return sqrt((w*EPSILON*template_funcs::DSQR(nc))/(betaH[i]*deff(i,t)));
				return sqrt( ( w * EPSILON * nc_sqr ) / ( betaH[i] * deff(i,t) ) );
			}
		}	
	}
	catch(std::range_error &e){
		std::cerr<<e.what();
		return 0; 
	}
}

double slab_wg::conf_fact(int i,bool t)
{
	// Confinement factor of mode in slab waveguide

	try{

		if(i<0 || i > nbeta(t) ){
			throw std::range_error("Error: double slab_wg::conf_fact(int i, bool t)\n Attempting to access arrays out of range\n"); 
			return 0;
		}
		else{
			if(t){ 
				double hh_sqr = template_funcs::DSQR(h(i,t)); 
				double pp = p(i, t); 
				double pp_sqr = template_funcs::DSQR(pp); 
				double qq = q(i, t);
				double qq_sqr = template_funcs::DSQR(qq); 
				
				//return (d+(1.0/(p(i,t)))*(1.0/( 1.0 + template_funcs::DSQR(h(i,t))/template_funcs::DSQR(p(i,t)) ) )+(1/(q(i,t)))*(1/(1+template_funcs::DSQR(h(i,t))/template_funcs::DSQR(q(i,t)))))/(deff(i,t));
				return ( d + (1.0/pp)*(1.0/( 1.0 + hh_sqr/pp_sqr ) ) + (1.0/qq)*(1.0/( 1.0 + hh_sqr/qq_sqr ) ) ) / ( deff(i,t) );
			}
			else{
				double hh_sqr = template_funcs::DSQR(h(i,t));
				double pp = p(i, t); 
				double pp_sqr = template_funcs::DSQR(pp); 
				double qq = q(i, t);
				double qq_sqr = template_funcs::DSQR(qq); 
				return ( d + (1.0/( g1(i,t) * pp ) )*(1.0/( 1.0 + hh_sqr/pp_sqr ) ) + ( 1.0 /( g2(i,t) * qq ) )*( 1.0 / (1.0+hh_sqr/qq_sqr ) ) ) / ( deff(i,t) );
			}
		}	
	}
	catch(std::range_error &e){
		std::cerr<<e.what();
		return 0; 
	}
}

double slab_wg::eigeneqn(double x,bool t)
{
	// Dispersion equation for slab waveguide written as a sum of sine and cosine functions
	// The version in eigeneqn_3, written as a sum over atan(*) is preferred for fast root finding

	try{

		if(x>=lower && x<=upper){
			/*double h=sqrt(template_funcs::DSQR(k)*template_funcs::DSQR(nc)-x*x);
			double p=sqrt(x*x-template_funcs::DSQR(k)*template_funcs::DSQR(ns));
			double q=sqrt(x*x-template_funcs::DSQR(k)*template_funcs::DSQR(ncl));*/

			double x_sqr = template_funcs::DSQR(x); 
			double h=sqrt(k_sqr_nc_sqr - x_sqr );
			double p=sqrt(x_sqr - k_sqr_ns_sqr);
			double q=sqrt(x_sqr - k_sqr_ncl_sqr);

			if(t){//TE modes
				return ( template_funcs::DSQR(h) - (p*q) ) * sin( ( d * h ) ) - h * ( p + q ) * cos( ( d * h ) );
			}
			else{//TM modes
				return ( template_funcs::DSQR(h) - aa * bb * ( p * q ) ) * sin( ( d * h ) ) - ( h * aa * p + h * bb * q ) * cos( ( d * h ) );
			}		
		}
		else{
			return 0; 
			throw std::range_error("Error: double slab_wg::eigeneqn(double x,bool t)\n Attempting to compute dispersion equation outside valid range\n"); 
		}	
	}
	catch(std::range_error &e){
		std::cerr<<e.what();
		return 0; 
	}
}

//double slab_wg::eigeneqn_2(double x,bool t,int mm)
//{
//	// Another version of the dispersion equation
//	// The roots of this equation are equal to the reduced effective index value
//	// This means that the root search is confined to the interval (0, 1)
//	if(t){//TE modes
//		return 2.0*V*sqrt((1-x))-mm*PI-atan(sqrt(x/(1-x)))-atan((x+g)/(1-x));
//	}
//	else{//TM modes
//		return 2.0*V*sqrt((1-x))-mm*PI-atan(sqrt(aa*(x/(1-x))))-atan(bb*((x+g)/(1-x)));
//	}
//}

double slab_wg::eigeneqn_3(double x,bool t,int mm)
{
	// Using this version of the dispersion equation means that you don't have to do a bracketing step
	// which you would have to do if you used eigeneqn

	try{
	
		if(k_sqr_nc_sqr > k_sqr_ns_sqr){

			double h,p,q,tmp;

			double x_sqr = template_funcs::DSQR(x);

			//tmp=template_funcs::DSQR(k)*template_funcs::DSQR(nc)-x*x;
			tmp=k_sqr_nc_sqr-x_sqr;
			h=(tmp>0?sqrt(tmp):0.0);
	
			//tmp=x*x-template_funcs::DSQR(k)*template_funcs::DSQR(ns);
			tmp=x_sqr-k_sqr_ns_sqr;
			p=(tmp>0?sqrt(tmp):0.0);
	
			//tmp=x*x-template_funcs::DSQR(k)*template_funcs::DSQR(ncl);
			tmp=x_sqr-k_sqr_ncl_sqr;
			q=(tmp>0?sqrt(tmp):0.0);

			if(t){//TE modes
				return d*h-mm*PI-atan(q/h)-atan(p/h);
			}
			else{//TM modes
				return d*h-mm*PI-atan(bb*(q/h))-atan(aa*(p/h));
			}			
		}
		else{
			std::string reason = "Error: double slab_wg::eigeneqn_3(double x,bool t,int mm)\n";
			reason += "Input parameters not correctly defined\n";
			reason += "k_{0}^{2} n_{c}^{2} = " + template_funcs::toString(k_sqr_nc_sqr) + ", k_{0}^{2} n_{s}^{2} = " + template_funcs::toString(k_sqr_ns_sqr) + "\n"; 
			return 0; 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double slab_wg::TE_TM(double x, int i, bool mode)
{
	// Function which defines the shape of the modes
	// This method uses the stored computed propagation constants
	// Assumes that the slab wg has core region defined on -t < x < t, t = d/2

	try{

		if(i<0 || i>nbeta(mode) ){
			throw std::range_error("Error: double slab_wg::TE_TM(double x,int i,bool mode)\n Attempting to access arrays out of range\n"); 
			return 0;
		}
		else{
			double t=d/2;
		
			if(x>t){
				return norm_const(i,mode)*cos(h(i,mode)*t-phase(i,mode))*exp(-q(i,mode)*(x-t)); // Cladding
			}
			else if(x<-t){
				return norm_const(i,mode)*cos(h(i,mode)*t+phase(i,mode))*exp(p(i,mode)*(x+t)); // Substrate
			}
			else{
				return norm_const(i,mode)*cos(h(i,mode)*x-phase(i,mode)); // Core
			}			
		}	
	}
	catch(std::range_error &e){
		std::cerr<<e.what();
		return 0; 
	}

	// Assumes that the slab wg has core region defined on -d < x < 0, otherwise the same
	/*if(x>0){
		return norm_const(i,mode)*exp(-q(i,mode)*x);
	}
	else if(x<-d){
		return norm_const(i,mode)*(cos(h(i,mode)*d)+(q(i,mode)/h(i,mode))*sin(h(i,mode)*d))*exp(p(i,mode)*(x+d));
	}
	else{
		return norm_const(i,mode)*(cos(h(i,mode)*x)-(q(i,mode)/h(i,mode))*sin(h(i,mode)*x));
	}*/
}

void slab_wg::calculate_all_modes(int n, double Lx, std::string &storage_directory)
{
	//Compute the modes and output the results to a txt file
	
	// compute effective indices and write solution statistics to a file
	calculate_all_neffs(storage_directory); 
	
	// output computed mode shapes to a file
	output_modes(TE, n, Lx, storage_directory);

	// output computed mode shapes to a file
	output_modes(TM, n, Lx, storage_directory);

	std::cout<<"Calculation Complete\n";
	std::cout<<"\n";
}

void slab_wg::calculate_neff(bool mode)
{
	// compute the effective indices for a particular polarisation
	// output the results
	// R. Sheehan 18 - 7 - 2014

	// This is really just a fancy wrapper for neff_search
	// It adds nothing except the printed output to the screen
	// R. Sheehan 31 - 8 - 2016

	neff_search(mode); 

	std::string pol = (mode ? "TE" : "TM" );

	std::cout<<"There are "<<static_cast<int>(nbeta(mode))<<" calculated "+pol+" modes\n";
	for(int i=0;i<static_cast<int>(nbeta(mode));i++){
		std::cout<<"beta["<<i+1<<"] = "<<std::setprecision(15)<<_beta(i,mode)<<"\n";
	}
	std::cout<<"\n";
}

void slab_wg::calculate_all_neffs(std::string &storage_directory)
{
	// compute the effective indices for both polarisations in a slab wg
	// R. Sheehan 18 - 7 - 2014

	calculate_neff(TE); 

	calculate_neff(TM); 

	output_all_stats(storage_directory); 
}

void slab_wg::neff_search(bool mode)
{
	//This version solves the standard slab dispersion equation with a superior root finder
	//sic. Brent's Method NRinC ch. 9
	//R. Sheehan 28 - 5 - 2010

	try{
		
		if(lower < upper){

			int m;
			double b;

			std::vector<double> vec;

			for(m=0; m < M; m++){
				b = zbrent(lower, upper, EPS, mode, m);

				if(b>lower && b<upper){
					vec.push_back(b);
				}

			}

			if(mode){
				betaE=vec;
			}
			else{
				betaH=vec;
			}

			vec.clear();		
		}
		else{
			std::string reason = "Error: void slab_wg::neff_search(bool mode)\n";
			reason += "Search range is not correctly defined\n";
			reason += "lower = " + template_funcs::toString(lower, 4) + ", upper = " + template_funcs::toString(upper, 4) + "\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double slab_wg::zbrent(double x1,double x2,double tol,bool t,int mm)
{
	//Find the root of the function between x1 and x2 using brent's method
	//The root is refined to +/- tol
	//Seemingly one of the best methods to use

	//This will be used to compute the roots of eigeneq_3
	//R. Sheehan 28 - 5 - 2010

	try{
	
		bool c1 = fabs(x1 - x2) > 1.0e-9 ? true : false; // cannot have x1 == x2
		bool c2 = tol > 1.0e-16 ? true : false; 
		bool c3 = mm < M ? true : false; 

		if(c1 && c2 && c3){

			int iter;
	
			static const int ITMAX=100;//Used in rtbis, rtflsp, rtsec, zriddr

			double a=std::min(x1, x2), b = std::max(x1, x2), c = std::max(x1, x2), d, e, min1, min2;
			double fc,p,q,r,s,tol1,xm;
			double fa = eigeneqn_3(a,t,mm), fb = eigeneqn_3(b,t,mm);

			if((fa>0.0 && fb>0.0) || (fa<0.0 && fb<0.0)){
				std::cerr<<"Root must be bracketed in zbrent\n";	
			}
			fc=fb;
			for(iter=1;iter<=ITMAX;iter++){
				if((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0)){
					c=a;
					fc=fa;
					e=d=b-a;
				}
				if(fabs(fc)<fabs(fb)){
					a=b;
					b=c;
					c=a;
					fa=fb;
					fb=fc;
					fc=fa;
				}
				tol1=2.0*EPS*fabs(b)+0.5*tol;
				xm=0.5*(c-b);
				if(fabs(xm)<=tol1 || fb==0.0) return b;
				/*if(fabs(xm)<=tol1 || fb==0.0){
					std::cout<<"Brent's Method converged in "<<iter<<" iterations\n";
					return b;
				}*/
				if(fabs(e)>=tol1 && fabs(fa)>fabs(fb)){
					s=fb/fa;
					if(a==c){
						p=2.0*xm*s;
						q=1.0-s;
					}
					else{
						q=fa/fc;
						r=fb/fc;
						p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
						q=(q-1.0)*(r-1.0)*(s-1.0);
					}
					if(p>0.0) q=-q;
					p=fabs(p);
					min1=3.0*xm*q-fabs(tol1*q);
					min2=fabs(e*q);
					if(2.0*p<std::min(min1,min2)){
						e=d;
						d=p/q;
					}
					else{
						d=xm;
						e=d;
					}
				}
				else{
					d=xm;
					e=d;
				}
				a=b;
				fa=fb;
				if(fabs(d)>tol1){
					b+=d;	
				}
				else{
					b+=template_funcs::SIGN(tol1,xm);	
				}
				fb=eigeneqn_3(b,t,mm);
			}
			std::cerr<<"Maximum number of iterations exceeded in zbrent\n";
			return 0.0;
		}
		else{
			std::string reason = "Error: double slab_wg::zbrent(double x1,double x2,double tol,bool t,int mm)\n"; 
			if(!c1) reason += "Cannot have x1 = x2\nx1 = " + template_funcs::toString(x1) + ", x2 = "+ template_funcs::toString(x1) + "\n";
			if(!c2) reason += "Desired tolerance is less than smallest allowable EPS\n";
			if(!c3) reason += "mm >= M\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output( e.what() ); 
		exit(EXIT_FAILURE); 
	}
}

void slab_wg::output_modes(bool mode, int N, double Lx, std::string &storage_directory)
{
	//This function will calculate the solutions corresponding to each mode and output them to a file

	double xi;

	std::string filename; 
	std::string pol = (mode ? "TE" : "TM" );

	std::ofstream write;

	if( nbeta( mode ) > 0 ){

		std::vector< std::vector< double > > mat;

		mat.resize( nbeta( mode ) + 2); // We want to output M solutions plus the corresponding positions
		
		double dx = Lx/( static_cast<double>(N-1) ); 

		xi = -0.5*Lx;
		
		mat[0].resize(N + 2);
		
		for(int j = 1; j <= N; j++ ){
			mat[0][j] = xi;
			xi += dx;
		}

		for(int i=1 ; i <= nbeta(mode) ; i++){
			mat[i].resize(N + 2);		
			for(int j=1;j<=N;j++){		
				mat[i][j] = TE_TM( mat[0][j], i - 1, mode );
			}
		}

		// Output all the modes to a single file
		filename = storage_directory + pol + "_Mode_Profiles.txt";
		write.open(filename.c_str(),std::ios_base::out|std::ios_base::trunc);

		if(write.is_open()){

			for(int i=1;i<=N;i++){
				for(int j=0;j<(nbeta(mode)+1);j++)
					if(j==nbeta(mode))
						write<<std::setprecision(20)<<mat[j][i];
					else
						write<<std::setprecision(20)<<mat[j][i]<<" , ";
				write<<"\n";
			}

			write.close();

		}

		// Output the modes to separate files
		/*for(int j = 0; j < ( nbeta(mode) + 1 ); j++){

				if(j==0){
					filename = storage_directory + pol + "_Positions" + dottxt;
				}
				else{
					filename = storage_directory + pol + "_Mode_" + template_funcs::toString(j) + dottxt;	
				}

				write.open(filename.c_str(),std::ios_base::out|std::ios_base::trunc);

				if(write.is_open()){
				
					for(int i=1; i<=N; i++){
						write<<std::setprecision(15)<<mat[j][i]<<"\n";
					}

					write.close(); 
				}
			}

			mat.clear();

		}*/

		// output the sine-cosine form of the dispersion equation

		double db = ( ( upper - lower ) / ( 100 - 1 ) );

		filename = storage_directory + pol + "_Dispersion_Eqn.txt";
		write.open( filename.c_str(), std::ios_base::out|std::ios_base::trunc);

		for(int i=1;i<99;i++){
			write<<lower+i*db<<" , "<<std::setprecision(20)<<eigeneqn(lower+i*db,mode)<<"\n";
		}

		write.close();	

		// Output each of the atan dispersion equations to their own file
		// Not sure if this is necessary
		//if(nbeta(mode) > 0){
		//	//char file[100];
		//	for(int m = 0; m < nbeta(mode); m++){
		//		//sprintf_s(file,100,"TE_Eigeneqn_%d.txt",m+1);
		//		filename = storage_directory + pol + "_Dispersion_Eqn_"+template_funcs::toString(m+1)+dottxt;
		//		write.open(filename.c_str(),std::ios_base::out|std::ios_base::trunc);
		//		for(int i=1;i<99;i++){
		//			write<<lower+i*db<<" , "<<std::setprecision(20)<<eigeneqn_3(lower+i*db,mode,m)<<"\n";
		//		}
		//		write.close();
		//	}

	}
}

void slab_wg::output_all_stats(std::string &storage_directory)
{
	//This function outputs all the stats associated with a particular calculation

	std::string file; 

	file = storage_directory + "Slab_WG_Stats.txt"; 

	std::ofstream stats;
	stats.open(file.c_str(),std::ios_base::out|std::ios_base::trunc);

	stats<<"Output File for the slab waveguide calculator\n";
	stats<<"\n";
	stats<<"Waveguide width = "<<d<<" microns\n";
	stats<<"Wavelength = "<<l<<" microns\n";
	stats<<"Wavenumber = "<<k<<" inverse microns\n";
	stats<<"\n";
	stats<<"Core Index = "<<nc<<"\n";
	stats<<"Substrate Index = "<<ns<<"\n";
	stats<<"Cladding Index = "<<ncl<<"\n";
	stats<<"\n";
	stats<<"Numerical Aperture = "<<na<<"\n";
	stats<<"V-Parameter = "<<V<<"\n";
	stats<<"Number of Modes = "<<M<<"\n";
	stats<<"\n";

	if(nbeta(TE) > 0){

		output_stats(TE, stats);

	}

	if(nbeta(TM) > 0){

		output_stats(TM, stats);

	}

	stats.close();
}

void slab_wg::output_stats(bool mode, std::ofstream &file_obj)
{
	// output the results for a particular polarisation
	// R. Sheehan 18 - 7 - 2014

	std::string pol = (mode ? "TE" : "TM" ); 

	if(file_obj.is_open()){
		
		file_obj<<pol<<" Modes\n";
		file_obj<<"There are "<<nbeta(mode)<<" calculated "+pol+" modes\n";

		for(int i=0;i<nbeta(mode);i++){
			file_obj<<pol+"_{"<<i+1<<"} = "<<std::setprecision(10)<<_beta( i, mode)<<" , n_eff_{"<<i+1<<"} = "<<std::setprecision(10)<<(_beta( i, mode)/k)<<"\n";
		}
		file_obj<<"\n";
		for(int i=0;i<nbeta(mode);i++){
			file_obj<<"Confinement Factor "<<i+1<<" = "<<std::setprecision(10)<<conf_fact(i,mode)<<"\n";
		}
		file_obj<<"\n";
		for(int i=0;i<nbeta(mode);i++){
			file_obj<<"Normalisation Constant "<<i+1<<" = "<<std::setprecision(10)<<norm_const(i,mode)<<"\n";
		}
		file_obj<<"\n";
		for(int i=0;i<nbeta(mode);i++){
			file_obj<<"Phase "<<i+1<<" = "<<std::setprecision(10)<<phase(i,mode)<<"\n";
		}
		file_obj<<"\n";
		for(int i=0;i<nbeta(mode);i++){
			file_obj<<"Effective Width "<<i+1<<" = "<<std::setprecision(10)<<deff(i,mode)<<"\n";
		}
		file_obj<<"\n";
		for(int i=0;i<nbeta(mode);i++){
			file_obj<<"h_"<<i+1<<" = "<<std::setprecision(10)<<h(i,mode)<<"\n";
		}
		file_obj<<"\n";
		for(int i=0;i<nbeta(mode);i++){
			file_obj<<"p_"<<i+1<<" = "<<std::setprecision(10)<<p(i,mode)<<"\n";
		}
		file_obj<<"\n";
		for(int i=0;i<nbeta(mode);i++){
			file_obj<<"q_"<<i+1<<" = "<<std::setprecision(10)<<q(i,mode)<<"\n";
		}
		file_obj<<"\n";
	}
}