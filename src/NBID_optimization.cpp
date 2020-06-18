/**
 *
 * This is the beginning of an attempt to optimize the NBID model using RcppNumerical.
 * This setup works in some cases but is VERY sensitive to initial parameter estimates.
 * In the cases where the estimates are "bad", the optimization spins off toward infinity.
 * Note: This does also happen in the R optimization with L-BFGS if bounds on parameters
 * aren't set, but assuming (loose) parameter bounds are set in R, that implementation
 * always behaves reasonably, converging arbitrarily towards the true parameter values as
 * sample size increases. This implementation needs debugging to figure out what's going
 * on when the optimization careens off toward infinity but right now I'm stumped.
 * (06/16/2020)
 *
 */

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppNumerical.h>

using namespace Rcpp;
using namespace Numer;

/* 
 * Used for sanity checking likelihood calculation
 * This code is essentially copied into the NBID class
 */
//' Calculate negative log likelihood
//' 
//' @param x parameter vector
//' @param counts count vector
//' @param groups per-sample group label vector (0 or 1)
//' @param size_factors per-sample size factor
//' @param null_model boolean indicating whether or not to use a group coefficient
//' @name negativeLL_NBID
//' @export
// [[Rcpp::export]]
double negativeLL_NBID(Eigen::VectorXd x, Eigen::VectorXd counts, Eigen::VectorXd groups, Eigen::VectorXd size_factors, bool null_model) {
   int n = counts.size();
   Eigen::VectorXd dispersion(n);
   for(int i = 0; i < n; i++) {
      dispersion[i] = x[groups[i]];
   }
   Eigen::VectorXd mu(n);
   if(null_model) {
      mu.fill(exp(x[2]));
   } else {
      for(int i = 0; i < n; i++) {
         mu[i] = exp(x[2] + groups[i]*x[3]);
      }
   }
   Eigen::VectorXd result(n);
   result.fill(0);
   for(int i = 0; i < n; i++) {
      double log_gamma = 0;
      if(counts[i] > 0) {
         for(int j = 1; j <= counts[i]; j++) {
            log_gamma += log(j);
         }
      }
      result[i] += - log_gamma - (pow(dispersion[i],-1)+counts[i])*log(1+dispersion[i]*size_factors[i]*mu[i]) + counts[i]*log(size_factors[i]*mu[i]) + counts[i]*log(dispersion[i]);
      if(counts[i] > 0) {
         for(int j = 0; j < counts[i]; j++) {
            result[i] += log(j + pow(dispersion[i], -1));
         }
      }
   }
   double nll = -result.sum();
   return(nll);
}

/*
 * Used for sanity checking gradient calculation
 * This code is essentially copied into the NBID class
 */
//' Calculate gradient
//' 
//' @param x parameter vector
//' @param counts count vector
//' @param groups per-sample group label vector (0 or 1)
//' @param size_factors per-sample size factor
//' @param null_model boolean indicating whether or not to use a group coefficient
//' @name gradient_NBID
//' @export
// [[Rcpp::export]]
Eigen::VectorXd gradient_NBID(Eigen::VectorXd x, Eigen::VectorXd counts, Eigen::VectorXd groups, Eigen::VectorXd size_factors, bool null_model) {
   int n = counts.size();
   Eigen::VectorXd dispersion(n);
   for(int i = 0; i < n; i++) {
      dispersion[i] = x[groups[i]];
   }
   Eigen::VectorXd mu(n);
   if(null_model) {
      mu.fill(exp(x[2]));
   } else {
      for(int i = 0; i < n; i++) {
         mu[i] = exp(x[2] + groups[i]*x[3]);
      }
   }
   Eigen::VectorXd grad(x.size());
   grad.fill(0);
   for(int i = 0; i < n; i++) {
      double a_i = 1 + dispersion[i]*size_factors[i]*mu[i];
      double a_inv = 1/dispersion[i];
      double a_inv_sq = pow(a_inv, 2);
      double ln_a = log(a_i);
      double mu_diff = counts[i] - size_factors[i]*mu[i];
      double b = 0;
      if(counts[i] > 0) {
         for(int j = 0; j < counts[i]; j++) {
            b = b + 1/(j + a_inv);
         }
      }
      // groups[i] is the index of x {0, 1} to update for group-specific dispersion
      grad[groups[i]] += a_inv_sq*(ln_a - b) + mu_diff/(dispersion[i]*a_i);
      grad[2] += mu_diff/a_i;
      if(!null_model) {
         grad[3] += (groups[i]*mu_diff)/a_i;
      }
   }
   grad = -1*grad;
   return(grad);
}

/*
 * Called by optimization
 */
class NBID: public MFuncGrad
{
   private:
      Eigen::VectorXd counts;
      Eigen::VectorXd groups;
      Eigen::VectorXd size_factors;
      bool null_model;
      int param_no;
   public:
      NBID(Eigen::VectorXd counts_, Eigen::VectorXd groups_, Eigen::VectorXd size_factors_, bool null_model_, int param_no_) :
      counts(counts_), groups(groups_), size_factors(size_factors_), null_model(null_model_), param_no(param_no_) {}
      
      // x are initial parameter values, grad are the gradient updates
      double f_grad(Constvec& x, Refvec grad)
      {
         int n = counts.size();
         Eigen::VectorXd dispersion(n);
         for(int i = 0; i < n; i++) {
            dispersion[i] = x[groups[i]];
         }
         Eigen::VectorXd mu(n);
         if(null_model) {
            mu.fill(exp(x[2]));
         } else {
            for(int i = 0; i < n; i++) {
               mu[i] = exp(x[2] + groups[i]*x[3]);
            }
         }
         Eigen::VectorXd result(n);
         result.fill(0);
         for(int i = 0; i < n; i++) {
            double log_gamma = 0;
            if(counts[i] > 0) {
               for(int j = 1; j <= counts[i]; j++) {
                  log_gamma += log(j);
               }
            }
            result[i] += - log_gamma - (pow(dispersion[i],-1)+counts[i])*log(1+dispersion[i]*size_factors[i]*mu[i]) + counts[i]*log(size_factors[i]*mu[i]) + counts[i]*log(dispersion[i]);
            if(counts[i] > 0) {
               for(int j = 0; j < counts[i]; j++) {
                  result[i] += log(j + pow(dispersion[i], -1));
               }
            }
         }
         double nll = -result.sum();

         // calculate gradient
         for(int i = 0; i < param_no; i++) {
            grad[i] = 0;
         }
         // print gradient updates; these start correct (same as R) but somehow updates send the optimization somewhere different
         // Rcout << "Gradient starting:" << std::endl << "\t" << grad[0] << std::endl << "\t" << grad[1] << std::endl << "\t" << grad[2] << std::endl;
         // if(param_no == 4) {
         //   Rcout << "\t" << grad[3] << std::endl;
         //}
         for(int i = 0; i < n; i++) {
            double a_i = 1 + dispersion[i]*size_factors[i]*mu[i];
            double a_inv = 1/dispersion[i];
            double a_inv_sq = pow(a_inv, 2);
            double ln_a = log(a_i);
            double mu_diff = counts[i] - size_factors[i]*mu[i];
            double b = 0;
            if(counts[i] > 0) {
               for(int j = 0; j < counts[i]; j++) {
                  b = b + 1/(j + a_inv);
               }
            }
            // groups[i] is the index of x {0, 1} to update for group-specific dispersion
            grad[groups[i]] += a_inv_sq*(ln_a - b) + mu_diff/(dispersion[i]*a_i);
            grad[2] += mu_diff/a_i;
            if(!null_model) {
               grad[3] += (groups[i]*mu_diff)/a_i;
            }
         }
         grad = -1*grad;
         // Rcout << "Gradient ending:" << std::endl << "\t" << grad[0] << std::endl << "\t" << grad[1] << std::endl << "\t" << grad[2] << std::endl;
         // if(param_no == 4) {
         //   Rcout << "\t" << grad[3] << std::endl;
         //}
         return(nll);         
      }
};

//' Optimize the NBID model
//' 
//' @param theta_init initial parameter vector
//' @param counts count vector
//' @param groups per-sample group label vector (0 or 1)
//' @param size_factors per-sample size factor
//' @param null_model boolean indicating whether or not to use a group coefficient
//' @param maxit maximum optimization steps
//' @name optim_NBID
//' @export
// [[Rcpp::export]]
Rcpp::List optimize_NBID(Eigen::VectorXd theta_init, Eigen::VectorXd counts, Eigen::VectorXd groups, Eigen::VectorXd size_factors, bool null_model, int maxit) {
   double fopt;
   NBID f(counts, groups, size_factors, null_model, theta_init.size());
   int res = optim_lbfgs(f, theta_init, fopt, maxit);
   return Rcpp::List::create(
      Rcpp::Named("xopt") = theta_init,
      Rcpp::Named("fopt") = fopt,
      Rcpp::Named("status") = res
   );   
}
