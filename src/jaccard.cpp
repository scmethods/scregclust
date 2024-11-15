#include <Rcpp/Lightest>

using length_type = std::vector<int>::size_type;

static length_type len_intersect(const std::vector<int>& x, const std::vector<int>& y) {
	length_type i = 0;
	length_type j = 0;
	length_type result = 0;

	while (i < x.size() && j < y.size()) {
		if (x[i] < y[j]) {
			i++;
		} else if (y[j] < x[i]) {
			j++;
		} else {
			result++;

			i++;
			j++;
		}
	}

	return result;
}

//' Perform the computations for thresholded Jaccard distance
//'
//' @details
//' This function is optimized for sparse matrices and computes the pairwise
//' Jaccard distances between the rows of the input matrix. Note that the
//' actual distance is not saved. Instead, a threshold (`eps`) is supplied
//' and an indicator matrix is returned, with a one indicating that the
//' distance is smaller than `eps` (equivalently, the Jaccard similarity
//' is larger than `1 - eps`).
//'
//' @param gs a list of integer vectors, one for each row, giving the column
//'           indices of the non-zero elements of the row or `NULL` if the
//'           whole row is empty.
//' @param eps an upper bound on the Jaccard distance (`1 - eps` becomes a
//'            lower bound on the Jaccard similarity)
//'
//' @return A list with row and column indices in the #row x #row indicator
//'         matrix specifying which rows in the original matrix had a distance
//'         of at most `eps`.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List jaccard_indicator_comp(Rcpp::List gs, double eps) {
	const auto n = static_cast<length_type>(gs.length());

	if (eps > 1.0 || eps < 0.0) {
		Rcpp::stop("0 <= eps <= 1 needs to hold");
	}

	std::vector<std::vector<int>> varr;
	varr.reserve(n);
	std::transform(gs.begin(), gs.end(), std::back_inserter(varr),
				   Rcpp::as<std::vector<int>>);

	const auto eps_ = 1.0 - eps;

	std::vector<int> ipairs;
	std::vector<int> jpairs;

	for (length_type i = 1; i < n; i++) {
		for (length_type j = 0; j < i; j++) {
			const auto len_inter = len_intersect(varr[i], varr[j]);
			const auto len_union = varr[i].size() + varr[j].size() - len_inter;

			if (static_cast<const double>(len_union) * eps_ <
				static_cast<const double>(len_inter)) {
				ipairs.emplace_back(i + 1);
				jpairs.emplace_back(j + 1);
			}
		}
	}

	Rcpp::List out;
	out["i"] = ipairs;
	out["j"] = jpairs;

	return out;
}
