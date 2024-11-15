#include <Rcpp/Lightest>

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include <algorithm>
#include <list>
#include <sstream>
#include <vector>

using Arr1d = Eigen::ArrayXd;
using Arr2d = Eigen::ArrayXXd;
using Matd = Eigen::MatrixXd;
using Vecd = Eigen::VectorXd;
using Veci = Eigen::VectorXi;

static Matd compute_xtx(const Matd& x) {
	const auto p = x.cols();

	Matd xtx = Eigen::MatrixXd::Zero(p, p);
	xtx.selfadjointView<Eigen::Lower>().rankUpdate(x.transpose());
	xtx.triangularView<Eigen::Upper>() = xtx.transpose();

	return xtx;
}

//' ADMM algorithm for solving the group-penalized least squares problem
//'
//' Implements estimation of the coop-lasso problem.
//'
//' @param y Target (n x m)
//' @param x Design matrix (n x p)
//' @param lambda Penalization parameter
//' @param weights A specific weight for each group (typically this is
//'                `sqrt(group size)`).
//' @param beta_0 Initial value for coefficients, allowing for warm start.
//'               Can be set to NULL, which results in the initial `beta`
//'               being a zero matrix.
//' @param rho_0 Initial ADMM step-size
//' @param alpha_0 Initial ADMM relaxation parameter
//' @param n_update Number of steps in-between updates of the
//'                 step-size/adaptation parameters
//' @param eps_corr Lower bound for the correlation in the step-size
//'                 update steps
//' @param max_iter Maximum number of iterations
//' @param eps_rel Relative tolerance for convergence check
//' @param eps_abs Absolute tolerance for convergence check
//' @param verbose Whether or not information about the optimization process
//'                should be printed to the terminal
//'
//' @return A list containing
//'     \item{beta}{The coefficients at convergence}
//'     \item{iterations}{Number of iterations}
//'
//' @references
//' Xu et al. (2017) Adaptive relaxed ADMM: Convergence theory and
//' practical implementation. DOI 10.1109/CVPR.2017.765
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List coop_lasso(
	Eigen::Map<Eigen::MatrixXd> y, Eigen::Map<Eigen::MatrixXd> x, double lambda,
	Eigen::Map<Eigen::ArrayXd> weights,
	Rcpp::Nullable<Rcpp::NumericMatrix> beta_0 = R_NilValue,  // Initialization
	double rho_0 = 0.2, double alpha_0 = 1.5, int n_update = 2,
	double eps_corr = 0.2,												 // Step-size
	int max_iter = 1000, double eps_rel = 1e-8, double eps_abs = 1e-12,	 // Convergence
	bool verbose = false) {
	// Record sizes
	const auto n = y.rows();
	const auto m = y.cols();
	const auto p = x.cols();

	if (x.rows() != n) {
		Rcpp::stop("y and x need to have the same number of rows.");
	}

	Matd beta = Eigen::MatrixXd::Zero(p, m);
	if (beta_0.isUsable()) {
		if (static_cast<Eigen::Index>(beta_0.as().nrow()) != p ||
			static_cast<Eigen::Index>(beta_0.as().ncol()) != m) {
			Rcpp::stop("beta_0 needs to be of size p x m");
		}

		beta =
			Eigen::Map<Matd>(beta_0.as().begin(), beta_0.as().nrow(), beta_0.as().ncol());
	}

	// Pre-compute some quantities to speed up computation
	// Precompute X^T X
	const Matd xtx = compute_xtx(x);
	// Precompute X^T Y
	const Matd xty = x.transpose() * y;

	// Ensure we are starting with a feasible point (i.e. zeta == beta)
	Matd zeta = beta;
	// For adaptive ADMM it is necessary to use unscaled multipliers
	// to get the computations right
	Matd mult = Eigen::MatrixXd::Zero(p, m);

	// We need to save one old iterate for step-size/relaxation
	// parameter estimation
	Matd beta_old = beta;
	Matd zeta_old = zeta;
	Matd mult_old = mult;
	Matd mult_hat_old = mult;

	// Initialize step-size and relaxation parameter
	if (rho_0 <= 0.0) {
		Rcpp::stop("rho_0 > 0 needs to hold");
	}
	auto rho = rho_0;
	double rho_old = 0.0;  // Ensure rho_old != rho in first iteration
	if (alpha_0 <= 0.0 || alpha_0 > 2.0) {
		Rcpp::stop("0 < alpha_0 < 2 needs to hold");
	}
	auto alpha = alpha_0;

	// Iteration counter
	int it = 1;
	// Parameter update counter
	int pc = 0;

	Matd xtx_rho_eye = xtx;
	Eigen::LDLT<Eigen::MatrixXd> ldlt;

	// ADMM algorithm
	while (true) {
		// Precompute if necessary
		if (pc == 0 || rho != rho_old) {
			// Only diagonal needs to be updated
			xtx_rho_eye.diagonal() = xtx.diagonal().array() + rho;
			ldlt.compute(xtx_rho_eye);
		}

		// Step 1: Update beta using pre-computed quantities
		beta = ldlt.solve(xty + rho * zeta + mult);

		// Relaxation step
		const Matd beta_relaxed = alpha * beta + (1.0 - alpha) * zeta;

		// Step 2: Update zeta
		Matd zeta_new = (beta_relaxed - mult / rho);

		const Arr1d shrink_pos =
			(1.0 -
			 lambda * weights / (rho * (zeta_new.cwiseMax(0.0).rowwise().norm()).array()))
				.max(0.0);
		const Arr1d shrink_neg =
			(1.0 -
			 lambda * weights / (rho * (zeta_new.cwiseMin(0.0).rowwise().norm()).array()))
				.max(0.0);

		for (Eigen::Index j = 0; j < m; j++) {
			for (Eigen::Index i = 0; i < p; i++) {
				if (zeta_new(i, j) >= 0.0) {
					zeta_new(i, j) *= shrink_pos(i);
				} else {
					zeta_new(i, j) *= shrink_neg(i);
				}
			}
		}

		// Step 3: Update multipliers
		mult += rho * (-beta_relaxed + zeta_new);

		// Convergence check
		// Compute primal and dual residuals
		const Matd primal_resid = -beta + zeta_new;
		const Matd dual_resid = rho * (zeta - zeta_new);

		if (verbose) {
			std::stringstream ss;
			ss << "\r#" << std::setw(5) << it << std::scientific << std::setprecision(4)
			   << " rho " << std::setw(5) << rho << " alpha " << std::setw(5) << alpha
			   << " prim_res " << std::setw(5) << primal_resid.norm() << " bnd "
			   << std::setw(5)
			   << fmax(eps_rel * fmax(beta.norm(), zeta_new.norm()), eps_abs)
			   << " dual_res " << std::setw(5) << dual_resid.norm() << " bnd "
			   << std::setw(5) << fmax(eps_rel * mult.norm(), eps_abs);

			Rcpp::Rcout << ss.str();
		}

		// Check residual convergence
		if ((primal_resid.norm() <=
			 fmax(eps_rel * fmax(beta.norm(), zeta_new.norm()), eps_abs)) &&
			(dual_resid.norm() <= fmax(eps_rel * mult.norm(), eps_abs))) {
			break;
		}

		// Step-size/relaxation parameter update
		pc++;
		if (pc == n_update) {
			// The hatted multipliers use non-relaxed beta and the zetas from
			// the previous iteration
			const Matd mult_hat = mult + rho * (-beta + zeta);

			const Matd delta_mult_hat = mult_hat - mult_hat_old;
			const Matd delta_h_hat = beta - beta_old;

			const Matd delta_mult = mult - mult_old;
			const Matd delta_g_hat = zeta_old - zeta;

			const auto norm_delta_mult_hat = delta_mult_hat.norm();
			const auto norm_delta_h_hat = delta_h_hat.norm();
			const auto norm_delta_mult = delta_mult.norm();
			const auto norm_delta_g_hat = delta_g_hat.norm();

			double a = 0.0;
			double a_corr = 0.0;

			if (norm_delta_mult_hat > 0.0 && norm_delta_h_hat > 0.0) {
				// Estimate local slope for h
				const auto delta_h_hat_delta_mult_hat =
					(delta_h_hat.array() * delta_mult_hat.array()).sum();
				const auto a_sd =
					delta_mult_hat.squaredNorm() / delta_h_hat_delta_mult_hat;
				const auto a_mg = delta_h_hat_delta_mult_hat / delta_h_hat.squaredNorm();

				if (2.0 * a_mg > a_sd) {
					a = a_mg;
				} else {
					a = a_sd - a_mg / 2.0;
				}

				a_corr =
					delta_h_hat_delta_mult_hat / (norm_delta_h_hat * norm_delta_mult_hat);
			}

			double b = 0.0;
			double b_corr = 0.0;

			if (norm_delta_mult > 0.0 && norm_delta_g_hat > 0.0) {
				// Estimate local slope for g
				const auto delta_g_hat_delta_mult =
					(delta_g_hat.array() * delta_mult.array()).sum();
				const auto b_sd = delta_mult.squaredNorm() / delta_g_hat_delta_mult;
				const auto b_mg = delta_g_hat_delta_mult / delta_g_hat.squaredNorm();

				if (2.0 * b_mg > b_sd) {
					b = b_mg;
				} else {
					b = b_sd - b_mg / 2.0;
				}

				b_corr = delta_g_hat_delta_mult / (norm_delta_g_hat * norm_delta_mult);
			}

			// Store old rho to check whether it changed and we need to
			// update pre-computed quantities
			rho_old = rho;

			// Update step-size if appropriate
			if (a_corr > eps_corr && b_corr > eps_corr) {
				rho = sqrt(a * b);
			} else if (a_corr > eps_corr && b_corr <= eps_corr) {
				rho = a;
			} else if (a_corr <= eps_corr && b_corr > eps_corr) {
				rho = b;
			}
			// Else: Leave rho as is

			// Update relaxation parameter if appropriate
			if (a_corr > eps_corr && b_corr > eps_corr) {
				alpha = 1.0 + 2.0 / (sqrt(a * b) * (1.0 / a + 1.0 / b));
			} else if (a_corr > eps_corr && b_corr <= eps_corr) {
				alpha = 1.9;
			} else if (a_corr <= eps_corr && b_corr > eps_corr) {
				alpha = 1.1;
			} else {
				alpha = 1.5;
			}

			// House-keeping
			beta_old = beta;
			zeta_old = zeta_new;
			mult_old = mult;
			mult_hat_old = mult_hat;

			// Reset counter
			pc = 0;
		}

		// House-keeping
		zeta = zeta_new;

		// Check iteration limit
		it++;
		if (it > max_iter) {
			if (verbose) {
				Rcpp::Rcout << std::endl;
			}
			Rcpp::Rcout << "Coop-Lasso: Maximum number of iterations reached";
			if (!verbose) {
				Rcpp::Rcout << std::endl;
			}
			break;
		}

		Rcpp::checkUserInterrupt();
	}

	if (verbose) {
		Rcpp::Rcout << std::endl;
	}

	Rcpp::List out;
	Rcpp::NumericMatrix beta_(beta.rows(), beta.cols(), beta.data());
	out["beta"] = beta_;
	out["iterations"] = it;

	return out;
}

// static void remove_kkt_elements(const Matd& beta, const Matd& grad, Matd& grad_bar) {
// 	const auto n = beta.rows();
// 	const auto m = beta.cols();

// 	for (Eigen::Index j = 0; j < m; j++) {
// 		for (Eigen::Index i = 0; i < n; i++) {
// 	 		if ((beta(i, j) == 0.0) && (grad(i, j) > 0.0)) {
// 	 			grad_bar(i, j) = 0.0;
// 	 		}
// 		}
// 	}
// }

static void greedy_coord_descent(const Matd& Q, Matd& beta, Matd& grad) {
	const auto n = beta.rows();
	const auto m = beta.cols();

	for (Eigen::Index t = 0; t < n; t++) {
		Eigen::Index empty_passive_sets = 0;

		for (Eigen::Index j = 0; j < m; j++) {
			// Determine maximum absolute gradient over passive set
			Eigen::Index p = -1;
			auto max_val = (((beta.col(j).array() > 0.0) || (grad.col(j).array() < 0.0))
								.cast<double>() *
							grad.col(j).array().abs())
							   .maxCoeff(&p);

			// Eigen::Index p = -1;
			// double max_val = 0.0;
			// for (Eigen::Index i = 0; i < n; i++) {
			// 	if ((beta(i, j) > 0.0) || (grad(i, j) < 0.0)) {
			// 		auto abs_grad = fabs(grad(i, j));
			// 		if (abs_grad > max_val) {
			// 			max_val = abs_grad;
			// 			p = i;
			// 		}
			// 	}
			// }

			// Perform coordinate descent on the selected coefficient
			if (max_val == 0.0) {
				empty_passive_sets++;
				continue;
			}

			const auto dbeta = fmax(0.0, beta(p, j) - grad(p, j) / Q(p, p)) - beta(p, j);
			beta(p, j) += dbeta;
			grad.col(j) += dbeta * Q.col(p);
		}

		if (empty_passive_sets == m) {
			break;
		}
	}
}

//' Compute NNLS coefficients
//'
//' Computes non-negative least squares coefficients with a matrix
//' right hand side.
//'
//' @param x Coefficient matrix (p x n matrix)
//' @param y Right hand side (p x m matrix)
//' @param eps Convergence tolerance
//' @param max_iter Maximum number of iterations
//'
//' @return A list containing
//'  \item{beta}{The estimated coefficient matrix}
//'  \item{iterations}{A vector containing the number of iterations needed
//'                    for the `i`-th column in `y` in the `i`-th entry.}
//'
//' @references
//' Duy Khuong Nguyen and Tu Bao Ho. Accelerated anti-lopsided algorithm
//' for nonnegative least squares. International Journal of Data Science
//' and Analytics, 3(1):23–34, 2017.
//'
//' Adapted from <https://github.com/khuongnd/nnls_antilopsided>
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List coef_nnls(Eigen::Map<Eigen::MatrixXd> x, Eigen::Map<Eigen::MatrixXd> y,
					 double eps = 1e-12, int max_iter = 1000L) {
	const auto n = x.cols();
	auto m = y.cols();	// Will be reduced whenever right-hand sides reach convergence

	// Pre-compute some quantities to speed up computation
	// Precompute X^T X
	const Matd xtx = compute_xtx(x);

	const Vecd inv_sqrt_diag_xtx = 1.0 / xtx.diagonal().array().sqrt();

	const Matd Q =
		xtx.array() * (inv_sqrt_diag_xtx * inv_sqrt_diag_xtx.transpose()).array();
	// Multiply -x^T y row-wise by the elements in inv_sqrt_diag_xtx
	Matd grad = (-x.transpose() * y).array().colwise() * inv_sqrt_diag_xtx.array();

	Matd beta_final = Eigen::MatrixXd::Zero(n, m);
	Matd beta = Eigen::MatrixXd::Zero(n, m);
	Matd grad_bar = grad;
	// remove_kkt_elements(beta, grad, grad_bar);
	grad_bar.array() *=
		(1 - ((beta.array() == 0.0) && (grad.array() > 0.0))).cast<double>();

	// Save necessary number of iterations
	std::vector<int> iterations(static_cast<std::vector<int>::size_type>(m));
	std::fill(iterations.begin(), iterations.end(), max_iter);

	std::list<int> remaining_obs(static_cast<std::list<int>::size_type>(m));
	std::iota(remaining_obs.begin(), remaining_obs.end(), 0);

	for (int l = 0; l < max_iter; l++) {
		const Matd beta_save = beta;
		const Matd grad_save = grad;

		// Exact line search algorithm over passive variables
		const Matd Q_grad_bar = Q * grad_bar;
		const Arr1d alpha1 = (grad_bar.colwise().squaredNorm()).array() /
							 (grad_bar.array() * Q_grad_bar.array()).colwise().sum();

		for (Eigen::Index j = 0; j < m; j++) {
			const auto a = alpha1(j);
			if ((a == a) && fabs(a) >= 1e-20 && fabs(a) < 1e30) {
				beta.col(j) -= a * grad_bar.col(j);
				grad.col(j) -= a * Q_grad_bar.col(j);
				for (Eigen::Index i = 0; i < n; i++) {
					if (beta(i, j) < 0.0) {
						// Correct for negative elements
						grad.col(j) -= beta(i, j) * Q.col(i);
						beta(i, j) = 0.0;  // Remove them from updated iterate
					}
				}
			}
		}

		// Greedy coordinate descent algorithm (First time)
		greedy_coord_descent(Q, beta, grad);

		// Accelerated search
		const Matd dbeta = beta_save - beta;
		const Matd Q_dbeta = Q * dbeta;
		const Arr1d alpha2 = (grad.array() * dbeta.array()).square().colwise().sum() /
							 (dbeta.array() * Q_dbeta.array()).colwise().sum();

		for (Eigen::Index j = 0; j < m; j++) {
			const auto a = alpha2(j);
			if ((a == a) && fabs(a) >= 1e-20 && fabs(a) < 1e30) {
				beta.col(j) -= a * dbeta.col(j);
				grad.col(j) -= a * Q_dbeta.col(j);
				for (Eigen::Index i = 0; i < n; i++) {
					if (beta(i, j) < 0) {
						// Correct for negative elements
						grad.col(j) -= beta(i, j) * Q.col(i);
						beta(i, j) = 0.0;  // Remove them from updated iterate
					}
				}
			}
		}

		// Greedy coordinate descent algorithm (Second time)
		greedy_coord_descent(Q, beta, grad);

		// Compute error
		grad_bar = grad;
		// remove_kkt_elements(beta, grad, grad_bar);
		grad_bar.array() *=
			(1 - ((beta.array() == 0.0) && (grad.array() > 0.0))).cast<double>();

		// Check for which rhs convergence has been achieved
		const Arr1d grad_bar_norms = grad_bar.colwise().norm();
		std::vector<Eigen::Index> kept_cols;
		kept_cols.reserve(remaining_obs.size());

		auto it = remaining_obs.begin();
		for (Eigen::Index i = 0; i < m; i++) {
			if (grad_bar_norms(i) < eps) {
				beta_final.col(*it) = beta.col(i);
				iterations[static_cast<std::vector<int>::size_type>(*it)] = l + 1;
				it = remaining_obs.erase(it);
			} else {
				kept_cols.push_back(i);
				it++;
			}
		}

		// Reduce problem to those rhs where convergence has not yet occurred
		m = static_cast<Eigen::Index>(kept_cols.size());
		if (m > static_cast<Eigen::Index>(0)) {
			// Use that columns in kept_cols are sorted by construction
			Eigen::Index j = 0;
			for (auto& i : kept_cols) {
				if (j != i) {
					beta.col(j) = beta.col(i);
					grad.col(j) = grad.col(i);
					grad_bar.col(j) = grad_bar.col(i);
				}
				j++;
			}

			beta.conservativeResize(Eigen::NoChange, m);
			grad.conservativeResize(Eigen::NoChange, m);
			grad_bar.conservativeResize(Eigen::NoChange, m);
		} else {
			break;
		}

		if (l == max_iter - 1) {
			Rcpp::Rcout << "NNLS: Maximum number of iterations reached" << std::endl;
		}

		Rcpp::checkUserInterrupt();
	}

	// Re-scale to original scale
	beta_final.array().colwise() *= inv_sqrt_diag_xtx.array();

	Rcpp::List out;
	Rcpp::NumericMatrix beta_final_(beta_final.rows(), beta_final.cols(),
									beta_final.data());
	out["beta"] = beta_final_;
	out["iterations"] = iterations;

	return out;
}
