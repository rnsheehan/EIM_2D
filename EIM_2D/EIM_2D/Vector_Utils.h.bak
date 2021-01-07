#ifndef VECTOR_UTILS_H
#define VECTOR_UTILS_H

// Declaration of a name space in which several useful std::vector related functions are defined
// R. Sheehan 7 - 12 - 2018

namespace vecut {
	std::vector<double> get_row(std::vector<std::vector<double>> &data, int r_num);
	std::vector<double> get_col(std::vector<std::vector<double>> &data, int c_num);

	void read_into_vector(std::string &filename, std::vector<double> &data, int &n_pts, bool loud = false);
	void write_into_file(std::string &filename, std::vector<double> &data, bool loud = false); 
	void read_into_matrix(std::string &filename, std::vector<std::vector<double>> &data, int &n_rows, int &n_cols, bool loud = false);

	void print_mat(std::vector<std::vector<double>> &matrix);
	void print_cmat(std::vector<std::vector<std::complex<double>>> &matrix);

	std::vector<std::vector<double>> mat_mat_product(std::vector<std::vector<double>> &mat1, std::vector<std::vector<double>> &mat2);
	std::vector<std::vector<std::complex<double>>> cmat_cmat_product(std::vector<std::vector<std::complex<double>>> &mat1, std::vector<std::vector<std::complex<double>>> &mat2);
	std::vector<std::complex<double>> cmat_cvec_product(std::vector<std::vector<std::complex<double>>> &mat, std::vector<std::complex<double>> &vec);

	std::vector<std::vector<double>> zero_mat(int& rows, int& cols);

	std::vector<std::vector<std::complex<double>>> zero_cmat(int &rows, int &cols); 
	std::vector<std::vector<std::complex<double>>> idn_cmat(int &size);
}

#endif
