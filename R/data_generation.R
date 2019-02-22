#' Simulate MAP data set:
#' @param G Number of genes to simulate
#' @param n_control Number of time points in control group
#' @param n_treat Number of time points in treatment group
#' @param n_rep Number of replicates in both group
#' @param k_real Always be 4
#' @param sigma2_r Variance parameter for \eqn{\tau_{g}}
#' @param sigma1_2_r Variance parameter for \eqn{\eta_{g1}}
#' @param sigma2_2_r Variance parameter for \eqn{\eta_{g2}}
#' @param mu1_r Mean parameter for \eqn{\eta_{g1}}
#' @param phi_g_r Dispersion parameter
#' @param p_k_real True proportion for each mixture component
#' @param x Time structured design for the simulated data
#' @details The vector of read counts for gene g, treatment group i, replicate j,
#' at time point \eqn{t,Y_{gij}(t)}, follows a Negative Binomial distribution
#' parameterized mean lambda_{gi} and phi_g, where
#' \eqn{E[Y_{gij}(t)] = lambda_{gi}(t)}.
#' \eqn{\lambda_{gi}(t)} is further modeled as
#' \eqn{\lambda_{gi}(t) = S_{ij} \exp[\eta_{g1}I_{i = 2} + B'(t)\eta_{g2}I_{i = 2} + B'(t)\tau_{g}]}
#' We have B'(t) are design matrix, which is constructed by 2 orthorgonal polynomial bases.
#' * t = 1,..., n_treat (or n_control if control group);
#' * j = 1,..., n_rep;
#' * g = 1,...,G; and
#' * \eqn{[\eta_{g1}, \eta_{g2}, \tau_{g}]} ~ 4-component gausssian mixture model
#' @return * Y1 Simulated data
#' @examples rm(list = ls())
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
#'                      phi_g_r = rep(1, 100),
#'                      p_k_real = c(0.7, 0.1, 0.1, 0.1),
#'                      x = x)
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
####################
data_generation = function(
  G = 100,
  n_control = 10,
  n_treat   = 10,
  n_rep     = 3,
  k_real = 4,
  sigma2_r = rep(1, 2),
  sigma1_2_r = 1,
  sigma2_2_r = c(3,2),
  mu1_r = 4,
  phi_g_r = rep(1, 100),
  p_k_real = c(0.7, 0.1, 0.1, 0.1),
  x = x #design matrix
){

  # dd: index of which mixture the data is simulated from
  dd = rep(c(0:(k_real-1)), p_k_real * G)

  ####basis function#####
  n_t <- n_rep * n_treat + n_rep * n_control
  X1 = c(rep(0, n_rep * n_control), rep(1, n_rep * n_treat))
  X2 = x
  n_basis = dim(X2)[2]
  eta_set = matrix(0, G, (2 * n_basis + 1))
  mu_set = matrix(0, G, n_t)
  Y1 = matrix(0, G, n_t)
  aa = rep(0, G)
  for (gg in 1:G){

    kk = dd[gg]+1
    eta2_coef = rnorm(n_basis) * sqrt(sigma2_r)

    if(kk %in% c(1,2)){
      eta1 = 0 #no scale difference

    }else{
      eta1 = rnorm(1, mu1_r, sd = sqrt(sigma1_2_r))

    }
    eta_set[gg,1] = eta1

    if(kk %in% c(1, 3)){
      eta3_coef = rep(0, n_basis) #not NPDE
    }else{
      eta3_coef = rnorm(n_basis) * sqrt(sigma2_2_r)

    }
    eta2 = x %*% eta2_coef
    eta_set[gg, 2 : (n_basis + 1)] = eta2_coef
    eta3 = X1 * x %*% eta3_coef
    eta_set[gg, (n_basis + 2): (2 * n_basis+1)] = eta3_coef

    mu_set[gg,] = exp(X1 * eta1 + eta2 + eta3)
    Y1[gg,] = rnbinom(n = n_t, mu = mu_set[gg,], size = 1/phi_g_r[gg])
  }

  return(Y1)


}



