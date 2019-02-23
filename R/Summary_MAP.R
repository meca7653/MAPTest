#' Summary for MAP Test
#' @param test_result MAP_test result:\cr
#'  A list of results for DE analysis
#' @param alpha Nominal level of FDR
#' @details The vector of read counts for gene g, treatment group i, replicate j,
#' at time point \eqn{t,Y_{gij}(t)}, follows a Negative Binomial distribution
#' parameterized mean \eqn{\lambda_{gi}} and \eqn{\phi_g}, where
#' \eqn{E[Y_{gij}(t)] = lambda_{gi}(t)}.
#' \eqn{\lambda_{gi}(t)} is further modeled as
#' \eqn{\lambda_{gi}(t) = S_{ij} exp[\eta_{g1}I_{i = 2} + B'(t)\eta_{g2}I_{i = 2} + B'(t)\tau_{g}]}
#' We have B'(t) are design matrix, which is constructed by 2 orthorgonal polynomial bases.
#' * t = 1,..., n_treat (or n_control if control group);
#' * j = 1,..., n_rep;
#' * g = 1,...,G; and
#' * \eqn{[\eta_{g1}, \eta_{g2}, \tau_{g}]} ~ 4-component gausssian mixture model.
#' After the parameter estimation we run the MAP test to detect the DE genes. At given
#' nominal level, a list of differential genes are returned at a given hypothesis.
#' @return A list of results for DE analysis.\cr
#' * Type Type of hypothesis
#' * Reject_index DE genes list
#' * FDR_hat Estimated FDR
#'
#' @examples
#' library(matlib)
#' set.seed(1)
#' n_basis = 2
#' n_control = 10
#' n_treat   = 10
#' n_rep = 3
#' tt_treat  = c(1:n_treat)/n_treat
#' nt = length(tt_treat)
#' ind_t = sort(sample(c(1:nt), n_control))
#' tt = tt_treat[ind_t]
#' tttt = c(rep(tt, n_rep), rep(tt_treat, n_rep))
#' z = x = matrix(0, length(tttt), n_basis)
#' z[,1] = 1.224745*tttt
#' z[,2] = -0.7905694 + 2.371708*tttt^2
#' x[,1] = z[,1] - Proj(z[,1], rep(1, length(tttt)))
#' x[,2] = z[,2] - Proj(z[,2], rep(1, length(tttt))) - Proj(z[,2], x[,1])
#' Y1 = data_generation(G = 100,
#'                      n_control = n_control,
#'                      n_treat   = n_treat,
#'                      n_rep     = n_rep,
#'                      k_real = 4,
#'                      sigma2_r = rep(1, 2),
#'                      sigma1_2_r = 1,
#'                      sigma2_2_r = c(3,2),
#'                      mu1_r = 4,
#'                      phi_g_r = rep(1,100),
#'                      p_k_real = c(0.7, 0.1, 0.1, 0.1),
#'                      x = x)
#' aaa <- proc.time()
#' est_result <- estimation_fun(n_control = n_control,
#'                              n_treat = n_treat,
#'                              n_rep = n_rep,
#'                              x = x,
#'                              Y1 = Y1,
#'                              nn = 300,
#'                              k = 4,
#'                              phi = NULL,
#'                              type = 2,
#'                              tttt = tttt)
#' aaa1 <- proc.time()
#' aaa1 - aaa
#'
#' G <- 100
#' k_real <- 4
#' p_k_real <- c(0.7, 0.1, 0.1, 0.1)
#' dd = rep(c(0:(k_real-1)), p_k_real * G)
#' result = MAP_test(est_result = est_result, Type = c(1:6), dd = dd, nn = 300)
#' ss <- Summary_MAP(result)
#' @import foreach
#' @import MASS
#' @import mvtnorm
#' @import randtoolbox
#' @import stats
#' @import mclust
#' @import EQL
#' @import matlib
#' @import parallel
#' @md
#' @export
Summary_MAP = function(test_result,
                    alpha = 0.05){
  result = list()
  Type = test_result$Type
  for(ii in c(1:length(Type))){
    FDR_hat <- test_result$FDR_final[max(which(test_result$FDR_final[,ii] < alpha)),ii]
    rej <- which(test_result$test_stat[,ii] < test_result$ct_final_all[
      max(which(test_result$FDR_final[,ii] < alpha)),ii])
    result[[ii]] = list(Type = Type[ii],
                      Reject_index = rej,
                      FDR_hat = FDR_hat)


  }
  return(result)
}





