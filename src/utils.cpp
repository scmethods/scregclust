#include <Rcpp/Lightest>

//' Allocate 3d-array and fill with matrix along first dimension
//'
//' @param input the matrix of size `n_obs x n_genes`
//' @param n_cl the size of the three-dimensional array's first dimension
//'
//' @return The allocated and filled array of size `n_cl x n_obs x n_genes`
//'
//' @keywords internal
// [[Rcpp::export]]
SEXP alloc_array(SEXP input, R_xlen_t n_cl) {
	const auto n_obs = static_cast<R_xlen_t>(Rf_nrows(input));
	const auto n_genes = static_cast<R_xlen_t>(Rf_ncols(input));

	const double* const pinput = REAL(input);

	const auto n_total = n_cl * n_obs * n_genes;
	if (n_total > R_XLEN_T_MAX) {
		Rcpp::stop("alloc_array: requested allocation too large");
	}

	SEXP arr = PROTECT(Rf_allocVector(REALSXP, n_total));
	double* const parr = REAL(arr);

	for (R_xlen_t i = 0, ub = n_obs * n_genes; i < ub; i++) {
		for (R_xlen_t j = 0; j < n_cl; j++) {
			parr[i * n_cl + j] = pinput[i];
		}
	}

	UNPROTECT(1);
	return arr;
}

//' Reset input 3d-array by filling matrix along first dimension
//'
//' @param arr The 3d-array of dimension `n_cl x n_obs x n_genes`
//' @param input The matrix of size `n_obs x n_genes`
//'
//' @keywords internal
// [[Rcpp::export]]
void reset_array(SEXP arr, SEXP input) {
	const int* const dims = INTEGER(PROTECT(Rf_getAttrib(arr, R_DimSymbol)));
	const auto n_cl = static_cast<R_xlen_t>(dims[0]);
	const auto n_obs = static_cast<R_xlen_t>(dims[1]);
	const auto n_genes = static_cast<R_xlen_t>(dims[2]);
	UNPROTECT(1);

	if (static_cast<R_xlen_t>(Rf_nrows(input)) != n_obs ||
		static_cast<R_xlen_t>(Rf_ncols(input)) != n_genes) {
		Rcpp::stop("reset_array: input has wrong dimensions");
	}

	const double* const pinput = REAL(input);
	double* const parr = REAL(arr);

	for (R_xlen_t i = 0, ub = n_obs * n_genes; i < ub; i++) {
		for (R_xlen_t j = 0; j < n_cl; j++) {
			parr[i * n_cl + j] = pinput[i];
		}
	}
}