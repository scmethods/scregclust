#include <Rcpp/Lighter>

#include <Eigen/Core>

// Assumes that
// - update_order is a permutation of 0:(length(k_) - 1)
// - indices in prior_indicators elements are within 0:(length(k_) - 1)
// [[Rcpp::export]]
Rcpp::IntegerVector allocate_clusters(SEXP resid_array,
									  Eigen::Map<Eigen::MatrixXd> resid_var,
									  Rcpp::List prior_indicator, Rcpp::IntegerVector k_,
									  Rcpp::IntegerVector update_order,
									  double prior_baseline, double prior_weight) {
	double* arr = REAL(resid_array);
	int* dims = INTEGER(PROTECT(Rf_getAttrib(resid_array, R_DimSymbol)));
	const auto n_cl = static_cast<Eigen::Index>(dims[0]);
	const auto n_obs = static_cast<Eigen::Index>(dims[1]);
	// const auto n_genes = static_cast<Eigen::Index>(dims[2]);
	UNPROTECT(1);

	const auto n_total = n_cl * n_obs;

	Rcpp::IntegerVector k(k_ - 1);

	Eigen::ArrayXd cluster_totals = Eigen::ArrayXd::Zero(n_cl);
	for (const auto& idx : k) {
		if (idx != -2) {
			cluster_totals[idx] += 1;
		}
	}

	// Iterate over genes in given order
	for (const auto& j : update_order) {
		// Load residuals for current gene
		const Eigen::Map<Eigen::MatrixXd> resid(arr + j * n_total, n_cl, n_obs);

		// Compute fraction of genes that gene j interacts with in each cluster,
		// according to prior information
		const Rcpp::IntegerVector prior_indices =
			Rcpp::as<Rcpp::IntegerVector>(prior_indicator[j]);
		Eigen::ArrayXd prior_frac = Eigen::ArrayXd::Zero(n_cl);
		for (const auto& idx : prior_indices) {
			if (k[idx] != -2) {
				prior_frac[k[idx]] += 1;
			}
		}

		for (Eigen::Index idx = 0; idx < n_cl; idx++) {
			if (cluster_totals[idx] > 0) {
				prior_frac[idx] /= cluster_totals[idx];
			}
		}
		// Add baseline to avoid numerical problems if a cluster is empty
		prior_frac += prior_baseline;
		// Convert to probabilities
		Eigen::ArrayXd prior_log_prob = prior_frac.log() - log(prior_frac.sum());

		// Compute model likelihood from residuals
		Eigen::MatrixXd model_log_likelihood =
			((-resid).array().colwise() / (2.0 * resid_var.row(j).transpose().array()))
				.colwise() -
			0.5 * (2.0 * M_PI * resid_var.row(j).transpose()).array().log();

		// Normalise to convert to probabilities
		Eigen::RowVectorXd max_model_ll = model_log_likelihood.colwise().maxCoeff();
		Eigen::MatrixXd model_ll_minus_max =
			model_log_likelihood.rowwise() - max_model_ll;
		Eigen::MatrixXd model_log_prob =
			model_ll_minus_max.rowwise() -
			model_ll_minus_max.array().exp().colwise().sum().log().matrix();

		// Compute total scores by weighting the model likelihood
		// and the prior probabilities. The latter are computed
		// by subtracting the total of the prior fractions.
		const Eigen::MatrixXd total_model_log_scores =
			(((1.0 - prior_weight) * model_log_prob.array()).colwise() +
			 (prior_weight * prior_log_prob));

		// Compute votes for each of the n_cl clusters
		Eigen::ArrayXi votes = Eigen::ArrayXi::Zero(n_cl);
		for (Eigen::Index i = 0; i < n_obs; i++) {
			int max_idx = -1;
			total_model_log_scores.col(i).maxCoeff(&max_idx);
			votes[max_idx] += 1;
		}

		// Move gene j to the new best cluster
		int best_cl = -1;
		votes.maxCoeff(&best_cl);
		if (k[j] != -2) {
			cluster_totals[k[j]] -= 1;
		}
		cluster_totals[best_cl] += 1;
		k[j] = best_cl;

		Rcpp::checkUserInterrupt();
	}

	return k + 1;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix compute_votes(SEXP resid_array, Eigen::Map<Eigen::MatrixXd> resid_var,
								  Eigen::Map<Eigen::MatrixXd> prior_model_log_prob) {
	double* arr = REAL(resid_array);
	int* dims = INTEGER(PROTECT(Rf_getAttrib(resid_array, R_DimSymbol)));
	const auto n_cl = static_cast<Eigen::Index>(dims[0]);
	const auto n_genes = static_cast<Eigen::Index>(dims[1]);
	const auto n_obs = static_cast<Eigen::Index>(dims[2]);
	UNPROTECT(1);

	const auto n_total = n_cl * n_genes;
	Rcpp::IntegerMatrix votes(static_cast<int>(n_genes), static_cast<int>(n_obs));

	for (Eigen::Index j = 0; j < n_obs; j++) {
		const Eigen::Map<Eigen::MatrixXd> resid(arr + j * n_total, n_cl, n_genes);

		const Eigen::MatrixXd model_log_likelihood =
			-0.5 * (2.0 * M_PI * resid_var.array()).log() -
			resid.transpose().array() / (2.0 * resid_var.array());
		const Eigen::MatrixXd total_model_log_scores =
			model_log_likelihood + prior_model_log_prob;

		for (Eigen::Index i = 0; i < n_genes; i++) {
			int max_idx;
			total_model_log_scores.row(i).maxCoeff(&max_idx);
			votes(i, j) = max_idx + 1;
		}

		Rcpp::checkUserInterrupt();
	}

	return votes;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix count_votes(Rcpp::IntegerMatrix votes, int n_cl) {
	Rcpp::IntegerMatrix counts(votes.nrow(), n_cl + 1);

	for (int j = 0; j < votes.ncol(); j++) {
		for (int i = 0; i < votes.nrow(); i++) {
			const auto v = votes(i, j);
			counts(i, v == -1 ? 0 : v) += 1;
		}
	}

	return counts;
}
