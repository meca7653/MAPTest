#' MAP test
#' @param est_result estimation_fun result:\cr
#' A List contains data_use (the data information) and result1 (Parameters estimations)
#' @param dd True index to check the true FDR
#' @param Type Type of hypothesis
#' * Type = 1 general DE detection
#' * Type = 2 PDE with no time-by-treatment and NPDE with both treatment and time-by-treatment
#' * Type = 3 NPDE with or without time-by-treatment
#' * Type = 4 PDE with no time-by-treatment
#' * Type = 5 NPDE with only time-by-treatment
#' * Type = 6 NPDE with both treatment and time-by-treatment
#' @param nn Number of QMC nodes used to estimate the test statistics
#' @details The vector of read counts for gene g, treatment group i, replicate j,
#' at time point \eqn{t,Y_{gij}(t)}, follows a Negative Binomial distribution
#' parameterized mean \eqn{\lambda_{gi}} and \eqn{\phi_g}, where
#' \eqn{E[Y_{gij}(t)] = \lambda_{gi}(t)}.
#' \eqn{\lambda_{gi}(t)} is further modeled as
#' \eqn{\lambda_{gi}(t) = S_{ij} \exp[\eta_{g1}I_{i = 2} + B'(t)\eta_{g2}I_{i = 2} + B'(t)\tau_{g}].}
#' We have \eqn{B'(t)} are design matrix, which is constructed by 2 orthorgonal polynomial bases.
#' * t = 1,..., n_treat (or n_control if control group);
#' * j = 1,..., n_rep;
#' * g = 1,...,G; and
#' * \eqn{[\eta_{g1}, \eta_{g2}, \tau_{g}]} ~ 4-component gausssian mixture model.
#' After the parameter estimation we run the MAP test to detect the DE genes.
#' @return A list of results for DE analysis.\cr
#' * Type Type of hypothesis of interest
#' * FDR_final FDR estimation of each hypothesis
#' * FDR_real The true FDR if dd is known
#' * power Power if dd is known
#' * T_x Likelihood result for each component
#' * test_stat Test statistics for each hypothesis
#' * ct_final_all Test statistics cut-off value at each estimated FDR level
#' @import foreach
#' @import MASS
#' @import mvtnorm
#' @import randtoolbox
#' @import stats
#' @import mclust
#' @import EQL
#' @import parallel
#' @import matlib
#' @examples
#' library(matlib)
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
#' @export
#' @md
#'
MAP_test = function(est_result,
                    dd = NULL,
                    Type = c(1:6),
                    nn = 6000){
  data_use <- est_result$data_use
  result1 <- est_result$result1
  n_basis <- data_use$n_basis
  # aa <- result1$aa
  # aa_use <- aa
  # aa_use[which(is.na(log(aa)))] = mean(aa, na.rm = T)
  phi <- 1/result1$phi
  k <- data_use$k
  X1 <- data_use$X1
  x <- data_use$x
  tttt <- data_use$tttt
  G <- dim(data_use$Y1)[1]
  Y1 <- data_use$Y1
  qmc.grid = halton(nn, 2 * n_basis+1, init = TRUE, normal = FALSE, usetime = FALSE)
  eta1_pre = qnorm(qmc.grid[,1])
  eta2_pre = t(qnorm(qmc.grid[,2:(n_basis + 1)]))
  eta3_pre = t(qnorm(qmc.grid[,(n_basis + 2) : (2 * n_basis + 1)]))
  wt = rep(1/nn,nn)
  eta2 = x %*% ((eta2_pre ) * sqrt(result1$sigma2))
  T_stat_fun = function(sigma2, sigma1_2 , sigma2_2, mu1 , p_k, phi){
    # phi = rep(1/aa_use, each = n_control + n_treat)
    # phi = 1/aa_use
    # phi <- 1/1

    mm_c_full = list()
    for (kk in c(1:k)){
      if(kk %in% c(1,2)){
        eta1 = rep(0, nn)

      }else{
        eta1 = eta1_pre * sqrt(sigma1_2) + mu1
      }

      if(kk %in% c(1,3)){
        eta3 = X1 * x %*% (eta3_pre) * 0
      }else{
        eta3 = X1 * x %*% (eta3_pre * sqrt(sigma2_2))
      }

      lambda = matrix(0, nrow = nn, ncol = length(tttt))

      for(ii in c(1:nn)){
        lambda[ii,] = exp(X1 * eta1[ii] + eta2[,ii] + eta3[,ii])
      }
      f_y_g_eta = function(ii, lambda){
        Y = c(t(Y1))
        lambda_use = c(replicate(G, lambda[ii,]))


        # res = apply(matrix(dnbinom(Y, size = phi, mu = lambda_use, log = T), ncol = G), 2, sum)
        res1 = apply(matrix((Y * (log(lambda_use) - log(phi + lambda_use)) +
                               phi * (log(phi) - log(phi + lambda_use))), ncol = G), 2, sum)
        return(res1)
      }


      res_pre <- foreach(ii = c(1:nn), .combine = "cbind")%do%{
        f_y_g_eta(ii, lambda = lambda)
      }
      res = apply(exp(res_pre), 1, mean)
      mm_c_full[[kk]] = res_pre
      # print(kk)

    }

    res = matrix(0, nrow = G, ncol = 4)
    res_full = matrix(0, nrow = G * length(tttt), ncol = 4)
    for(kk in c(1:4)){
      res[, kk] = log(apply(exp(mm_c_full[[kk]]), 1, mean)) + log(p_k[kk])
      res[res[, kk] == -Inf, kk] = -2000
    }
    return(res)
  }

  T_x = T_stat_fun(sigma2 = result1$sigma2,
                   sigma1_2 = result1$sigma1_2,
                   sigma2_2 = result1$sigma2_2,
                   mu1 =  result1$mu1,
                   p_k = result1$p_k,
                   phi = phi)


  T_x = (T_x - apply(T_x, 1, max))
  FDR_ct = seq(1e-5, 0.6, length.out = 100)
  if(is.null(dd)){
    FDR <- matrix(0, 100, length(Type))
    FDR_final = ct_final_all = matrix(0, length(FDR_ct), length(Type))
    test_stat = matrix(0, G, length(Type))
    for(ii in c(1:length(Type))){
      kk = Type[ii]
      # progress(ii/length(Type) * 100)
      # print(paste0(ii/length(Type) * 100, "%") )

      if(kk == 1){
        index_use = 1 #general DE detection
      }else if(kk == 2){
        index_use = c(1,2) #mean difference detection
      }else if(kk == 3){
        index_use = c(1,3) #NPDE detection
      }else if(kk == 4){
        index_use = c(1,2,4) #PDE with no time-by-treatment
      }else if(kk == 5){
        index_use = c(1,3,4) # NPDE with only time-by-treatment
      }else if(kk == 6){
        index_use = c(1,2,3)
      }
      # print(index_use)
      res = apply(t(t(exp(T_x)[,index_use])), 1, sum)/apply(exp(T_x), 1, sum)
      test_stat[,ii] = res
      ct =  seq(min(res),quantile(res)[4], length.out = 100)

      for (jj in c(1:length(ct))){
        phi_x = (res<=ct[jj])
        FDR[jj,ii] = sum(res*phi_x)/sum(phi_x)

      }
      for(jj in c(1:length(FDR_ct))){
        ct_final = max(ct[which(FDR[,ii] == max(FDR[which (FDR[,ii]<=FDR_ct[jj]), ii]))] )
        ct_final_all[jj, ii] = ct_final
        phi_x = (res<=ct_final)
        FDR_final[jj, ii] = sum(res*phi_x)/sum(phi_x)

      }


    }
    result = list(Type = Type,
                  FDR_final = FDR_final,
                  T_x = T_x,
                  test_stat = test_stat,
                  ct_final_all = ct_final_all)
  }else{
    FDR = matrix(0, 100, length(Type))
    FDR_final = FDR_real = power = ct_final_all = matrix(0, length(FDR_ct), length(Type))
    test_stat = matrix(0, G, length(Type))
    for(ii in c(1:length(Type))){
      # print(paste0(ii/length(Type) * 100, "%") )
      kk = Type[ii]
      # progress(ii/length(Type) * 100)
      if(kk == 1){
        index_use = 1 #general DE detection
      }else if(kk == 2){
        index_use = c(1,2) #mean difference detection
      }else if(kk == 3){
        index_use = c(1,3) #NPDE detection
      }else if(kk == 4){
        index_use = c(1,2,4) #PDE with no time-by-treatment
      }else if(kk == 5){
        index_use = c(1,3,4) # NPDE with only time-by-treatment
      }else if(kk == 6){
        index_use = c(1,2,3)
      }
      # print(index_use)
      res = apply(t(t(exp(T_x)[,index_use])), 1, sum)/apply(exp(T_x), 1, sum)
      test_stat[,ii] = res
      ct =  seq(min(res),quantile(res)[4], length.out = 100)

      for (jj in c(1:length(ct))){
        phi_x = (res<=ct[jj])
        FDR[jj,ii] = sum(res*phi_x)/sum(phi_x)

      }


      for(jj in c(1:length(FDR_ct))){
        ct_final = max(ct[which(FDR[,ii] == max(FDR[which (FDR[,ii]<FDR_ct[jj]), ii]))] )
        ct_final_all[jj, ii] = ct_final
        phi_x = (res<ct_final)
        FDR_final[jj, ii] = sum(res*phi_x)/sum(phi_x)
        FDR_real[jj, ii] = sum(phi_x[which(dd %in% c(index_use-1))])/sum(phi_x)
        power[jj, ii] = sum(phi_x[which(dd %in% c(0:(k-1))[-index_use])])/sum(dd %in%  c(0:(k-1))[-index_use])
      }
    }
    result = list(Type = Type,
                  FDR_final = FDR_final,
                  FDR_real = FDR_real,
                  power = power,
                  T_x = T_x,
                  test_stat = test_stat,
                  ct_final_all = ct_final_all)
  }
  return(result)
}





