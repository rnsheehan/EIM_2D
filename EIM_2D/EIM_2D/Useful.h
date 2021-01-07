#ifndef USEFUL_H
#define USEFUL_H

// Library of functions that are very useful
// R. Sheehan 4 - 7 - 2011

namespace useful_funcs{
	
	std::string TheTime();

	void exit_failure_output(std::string reason);

	//void read_into_vector(std::string &filename, std::vector<double> &data, int &n_pts, bool loud = false); 

	bool valid_filename_length(const std::string& name);
}

// Data type for an interval [xlower, xupper]
class interval {
public:
	// Constructor
	interval();
	interval(double xl, double xu);

	// Methods

	void set_xl_xu(double xl, double xu);

	bool has_bounds() { return interval_defined; }

	double get_x_lower() { return xlower; }
	double get_x_upper() { return xupper; }

private:
	bool interval_defined;
	double xlower;
	double xupper;
};


#endif