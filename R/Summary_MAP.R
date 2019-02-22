#' Summary for MAP Test
#' @param test_result MAP_test result:\cr
#'  A list of results for DE analysis
#' @param alpha Nominal level of FDR
#'
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
#'
#' @return A list of results for DE analysis.\cr
#' * Type Type of hypothesis
#' * Reject_index DE genes list
#' * FDR_hat Estimated FDR
#' @examples
#' aaa <- proc.time()
#' G <- 100
#' k_real <- 4
#' p_k_real <- c(0.7, 0.1, 0.1, 0.1)
#' dd = rep(c(0:(k_real-1)), p_k_real * G)
#' result <- MAP_test(est_result = est_result, Type = c(1:5), dd = dd, nn = 3000)
#' aaa1 <- proc.time()
#' aaa1 - aaa
#' Summary_MAP(result)
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





