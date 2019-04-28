#' Estimation for MAP data set:
#' @param n_control Number of time points in control group
#' @param n_treat Number of time points in treatment group
#' @param n_rep Number of replicates in both group
#' @param x Time structured design for the simulated data
#' @param Y1 RNA-seq data set: G by (n_control + n_treat) * n_rep matrix
#' @param k Always 4
#' @param nn Number of QMC nodes used in likelihood estimation
#' @param phi the dispersion parameter
#' @param type dispersion parameter estimation type
#' @param tttt Time points used in the design
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
#' We used latented negative binomial model with EM algorithm to estimate the paramters of mixture model.
#' @return A list of parameters that we need to use for DE analysis.\cr
#'
#' data_use List includes the data information
#' * Y1 RNA-seq data set: G by (n_control + n_treat) * n_rep matrix
#' * start Starting point of the EM algorithm
#' * k always 4
#' * n_basis Number of basis function have been used to construct the design matrix
#' * X1 Vector indicate control data (0) or treatment data (1)
#' * x Design matrix been used for estimation
#' * tttt Time points used in the design \cr
#'
#' result1 List includes the Parameters estimations
#' * mu1 Mean parameter for \eqn{\eta_{g1}}
#' * sigma1_2 Variance parameter for \eqn{\eta_{g1}}
#' * sigma2 Variance parameter for \eqn{\tau_{g}}
#' * sigma2_2 Variance parameter for \eqn{\eta_{g2}}
#' * p_k Proportion for each mixture component
#' * aa Dispertion parameter
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
#' #------------gaussian basis construction------------------
#' # method <- "gaussian"
#' # n_basis <- 2
#' # a = 1/4
#' # b = 2
#' # c = sqrt(a^2 + 2 * a * b)
#' # A = a + b + c
#' # B = b/A
#' # phi_cal = function(k, x){
#' #   Lambda_k =sqrt(2 * a / A) * B^k
#' #   1/(sqrt(sqrt(a/c) * 2^k * gamma(k+1))) *
#' #      exp(-(c-a) * x^2) * hermite(sqrt(2 * c) * x, k, prob = F)
#' #
#' # }
#' # z = do.call(cbind, lapply(c(1:n_basis), phi_cal, x = (tttt - mean(tttt))))
#' # x = matrix(0, length(tttt), n_basis)
#' # if(n_basis == 2){
#' #   x[,1] = z[,1] - Proj(z[,1], rep(1, length(tttt)))
#' #   x[,2] = z[,2] - Proj(z[,2], rep(1, length(tttt))) - Proj(z[,2], x[,1])
#' #
#' # }else{
#' #   x[,1] = z[,1] - Proj(z[,1], rep(1, length(tttt)))
#' #   x[,2] = z[,2] - Proj(z[,2], rep(1, length(tttt))) - Proj(z[,2], x[,1])
#' #   x[,3] = z[,3] - Proj(z[,3], rep(1, length(tttt))) - Proj(z[,3], x[,1]) - Proj(z[,3], x[,2])
#' # }
#' @import foreach
#' @import MASS
#' @import mvtnorm
#' @import randtoolbox
#' @import stats
#' @import mclust
#' @import EQL
#' @import parallel
#' @import matlib
#' @md
#' @export
estimation_fun = function(n_control = 10,
                          n_treat   = 10,
                          n_rep     = 3,
                          x,
                          Y1,
                          k = 4,
                          nn = 300,
                          phi = NULL,
                          type = NULL,
                          tttt
){
  aa_use = NULL
  gg = NULL
  n_basis <- dim(x)[2]
  x_full <- paste0(rep("X", 2 * n_basis + 1), c(1:(n_basis * 2 + 1)))
  x_reduce_all <- paste0(rep("X", n_basis), c(1:n_basis)+1)
  x_reduce <- paste0(rep("X", n_basis + 1), c(2:(n_basis*2 + 1)))
  x_reduce_eta_3 <- paste0(rep("X", n_basis + 1), c(1:(n_basis + 1)))
  FFormula <- as.formula(paste("Y1 ~", paste(x_full, collapse = "+"), "- 1"))
  FFormula_null <- as.formula(paste("Y1 ~", paste(x_reduce_all, collapse = "+"), "- 1"))
  FFormula_null_eta_1 <- as.formula(paste("Y1 ~", paste(x_reduce, collapse = "+"), "- 1"))
  FFormula_null_eta_3 <- as.formula(paste("Y1 ~", paste(x_reduce_eta_3, collapse = "+"), "- 1"))
  G <- dim(Y1)[1]
  X1 <- c(rep(0, n_rep * n_control), rep(1, n_rep * n_treat))
  #-----------GLM and M_cluster initialization----------------
  result <- foreach(gg = c(1:G), .combine = "rbind") %do% {
    print(gg)


    data = data.frame(cbind(Y1[gg, ], X1, x, X1 * x))
    names(data) = c("Y1", paste0(rep("X", n_basis + 1), c(1:(n_basis * 2 + 1))))


    possibleError <- tryCatch(
      glm.nb(FFormula, data = data, link = log),
      error = function(e)
        e
    )
    possibleError_1 <- tryCatch(
      glm.nb(FFormula_null, data = data, link = log),
      error = function(e)
        e
    )
    possibleError_2 <- tryCatch(
      glm.nb(FFormula_null_eta_3, data = data, link = log),
      error = function(e)
        e
    )


    ee_index = 0
    if(!inherits(possibleError, "error") | !inherits(possibleError_1, "error") | !inherits(possibleError_2, "error")){
      if (!inherits(possibleError, "error")) {
        #REAL WORK
        glm1 <- glm.nb(FFormula, data = data, link = log)

        coef <- glm1$coefficients
        vv_coef <- summary(glm1)$coef
        aa_glm <- glm1$theta


      }else{
        ee_index = 1
        # lm1 = lm(Y1 ~ . -1, data = data)
        coef = NA
        vv_coef = NA
        aa_glm <- NA

      }
      if(!inherits(possibleError_1, "error")) {

        glm2 <- glm.nb(FFormula_null, data = data, link = log)
        pp_v <- anova(glm1, glm2, test = "Chisq")$Pr[2]
        # print(pp_v)

      }else{
        pp_v <- NA
      }
      if(!inherits(possibleError_2, "error")) {
        glm3 <- glm.nb(FFormula_null_eta_3, data = data, link = log)
        pp_v_eta_3 <- anova(glm1, glm3, test = "Chisq")$Pr[2]

      } else{
        pp_v_eta_3 <- NA
      }
    }else{
      ee_index = 1
      # lm1 = lm(Y1 ~ . -1, data = data)
      coef = NA
      vv_coef = NA
      pp_v = NA
      pp_v_eta_3 = NA
      aa_glm <- NA

    }
    Y1_n <- Y1[gg,]
    m1_c <-
      apply(matrix(Y1_n[1:(n_control * n_rep)], nrow = n_control), 1, mean)
    m1_t <-
      apply(matrix(Y1_n[(1:(n_treat * n_rep)) + n_control * n_rep], nrow = n_treat), 1, mean)
    v1_c <-
      apply(matrix(Y1_n[1:(n_control * n_rep)], nrow = n_control), 1, var)
    v1_t <-
      apply(matrix(Y1_n[(1:(n_treat * n_rep)) + n_control * n_rep], nrow = n_treat), 1, var)
    aa <-
      sum(sum(v1_c - m1_c) + sum(v1_t - m1_t)) / (sum(m1_c ^ 2) + sum(m1_t ^ 2))


    list(coef, ee_index, vv_coef, pp_v, pp_v_eta_3, aa, aa_glm)

  }

  ee_index <- do.call(rbind, result[,2])
  coef <- do.call(rbind, result[,1])
  index_use <- which(ee_index == 0 & abs(coef[,1]) < 20)
  result_pre <- result
  result <- result_pre[index_use, ]
  vv_coef <- do.call(rbind, result[,3])
  aa <- do.call(rbind, result[,6])
  if(is.null(phi)){
    if(type == 1){
      aa_use <- aa
      aa_use[which(aa <= 0)] = mean(aa, na.rm = T)
      phi = rep(aa_use, each = n_control + n_treat)
    }else{
      aa_use <- mean(aa, na.rm = T)
      phi = aa_use
    }

  }else{
    phi = rep(1/phi, each = n_control + n_treat)
  }



  aa_glm <- do.call(rbind, result[,7])
  coef <- matrix(c(vv_coef[, 1]), nrow = 2 * n_basis + 1)
  eta_nb_use = t(coef)
  mmm1 = Mclust(eta_nb_use[,1], G = 2)
  mu1 = mmm1$parameters$mean[2]
  sigma1_2 = max(mmm1$parameters$variance$sigmasq)
  pp_v_eta_3 <- do.call(rbind, result[,5])
  QQ_2 = apply(eta_nb_use, 2, quantile)
  sigma2 = rep(0, n_basis)
  for(kk in c(1:n_basis) + 1){
    Low_2 = QQ_2[2,kk] - 1.5 * IQR(eta_nb_use[,kk])
    Upp_2 = QQ_2[4,kk] + 1.5 * IQR(eta_nb_use[,kk])
    ii_index = (eta_nb_use[,kk])[eta_nb_use[,kk] <= Upp_2 & eta_nb_use[,kk] >= Low_2]
    vvv = matrix(vv_coef[,1], nrow = 2 * n_basis + 1)
    sigma2[kk-1] = var(ii_index) +
      median((vvv[kk, ])[eta_nb_use[,kk] <= Upp_2 & eta_nb_use[,kk] >= Low_2 ]^2)

  }

  sigma2_2 = rep(0, n_basis)
  for(kk in c(1:n_basis) + 1 + n_basis){
    Low_2 = QQ_2[2,kk] - 1.5 * IQR(eta_nb_use[,kk])
    Upp_2 = QQ_2[4,kk] + 1.5 * IQR(eta_nb_use[,kk])
    ii_index = (eta_nb_use[,kk])[eta_nb_use[,kk] <= Upp_2 & eta_nb_use[,kk] >= Low_2]
    vvv = matrix(vv_coef[,1], nrow = 2 * n_basis + 1)
    sigma2_2[kk-1-n_basis] = var(ii_index) +
      median((vvv[kk, ])[eta_nb_use[,kk] <= Upp_2 & eta_nb_use[,kk] >= Low_2 ]^2)

  }
  p_k = rep(0, 4)
  p_k[1] = sum((pp_v_eta_3[mmm1$classification == 1])>0.1)/G
  p_k[2] = sum((pp_v_eta_3[mmm1$classification == 1])<=0.1)/G
  p_k[3] = sum((pp_v_eta_3[mmm1$classification == 2])>0.1)/G
  p_k[4] = sum((pp_v_eta_3[mmm1$classification == 2])<=0.1)/G
  if(sum(p_k, na.rm = T) < 1){
    p_k[which(is.na(p_k))] = (1 - sum(p_k, na.rm = T))/sum(is.na(p_k))
  }

  print(p_k)
  ######################################

  res_get = list(aa = aa, mu1 = mu1, sigma1_2 = sigma1_2,
                 sigma2 = sigma2, sigma2_2 = sigma2_2, p_k = p_k, phi = phi)
  data_use = list(Y1 = Y1, start = res_get, k = k,
                  n_basis = n_basis,
                  X1 = X1, x = x, tttt = tttt)
  #---------------------------Fitting-----------------------------------
  # aa_use = mean(aa, na.rm = T)


  #--------------------------Likelihood Function-------------------------
  likelihood_SP = function(sigma2, sigma1_2, mu1, p_k, sigma2_2, index, wt,
                           Y1 = Y1,
                           eta1_pre = eta1_pre,
                           eta2_pre = eta2_pre,
                           eta3_pre = eta3_pre,
                           k = k, X1 = X1, x = x, phi = 1/phi,
                           G = G, tttt = tttt){
    if(is.null(index)){
      index = list(c(1:nn), c(1:nn), c(1:nn), c(1:nn))
    }
    if(is.null(wt)){
      wt = list(rep(1/nn, nn), rep(1/nn, nn), rep(1/nn, nn), rep(1/nn, nn))
    }
    tmp = matrix(0, nrow = G, ncol = k)
    tmp_res = list()
    f_y_g_eta = function(ii, lambda){
      Y = c(t(Y1))
      lambda_use = c(replicate(G, lambda[ii,]))
      res1 = apply(matrix((Y * (log(lambda_use) - log(phi + lambda_use)) +
                             phi * (log(phi) - log(phi + lambda_use))), ncol = G), 2, sum)
      return(res1)
    }
    mc_y_g_eta = function(index, lambda, wt){

      res_pre = foreach(ii = index,
                        .combine = "cbind") %dopar%{
                          f_y_g_eta(ii, lambda)
                        }
      res = exp(res_pre)%*%wt
      return(list(res = log(res), res_pre = res_pre))
    }
    eta2 = x %*% ((eta2_pre ) * sqrt(sigma2))
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

      for(ii in index[[kk]]){
        lambda[ii,] = exp(X1 * eta1[ii] + eta2[,ii] + eta3[,ii])
      }

      mm_c = mc_y_g_eta(index = index[[kk]], lambda = lambda, wt = wt[[kk]])
      tmp[,kk] =  log(p_k[kk]) + mm_c$res
      tmp_res[[kk]] = mm_c$res_pre
      tmp[(tmp[,kk] == -Inf),kk] = -2000

    }

    ll = sum(tmp)

    return(list(ll = -(ll), res_pre = tmp_res, tmp = tmp))
  }

  likelihood_f_sig1 = function(par){
    sigma1_2 = par[1]
    mu1 = par[2]
    sigma2 = par[3:(n_basis+2)]
    sigma2_2 = par[(n_basis + 3):(2 * n_basis + 2)]
    print(par)
    tmp = likelihood_SP(sigma2, sigma1_2, mu1, p_k, sigma2_2, index, wt,
                        Y1 = Y1,
                        eta1_pre = eta1_pre,
                        eta2_pre = eta2_pre,
                        eta3_pre = eta3_pre,
                        k = k, X1 = X1, x = x, phi = 1/phi,
                        G = G, tttt = tttt)$tmp
    ll = sum(c(G_ugk * tmp)[is.na(c(G_ugk * tmp)) == FALSE])

    return(-(ll))
  }




  #-----------------------------QMC nodes-----------------------------
  qmc.grid = halton(nn, 2 * n_basis+1, init = TRUE, normal = FALSE, usetime = FALSE)
  dim(qmc.grid)

  eta1_pre = qnorm(qmc.grid[,1])
  eta2_pre = t(qnorm(qmc.grid[,2:(n_basis + 1)]))
  eta3_pre = t(qnorm(qmc.grid[,(n_basis + 2) : (2 * n_basis + 1)]))
  wt = rep(1/nn,nn)

  #----------------------------EM algorithm-------------------------
  index = NULL
  wt = NULL
  j = 0
  ll0 = 0
  ll1 = 100
  while (abs(ll0-ll1)/abs(ll0) >= 1e-2){
    j = j+1
    print(paste0("j_", j))
    ll0 = ll1
    p_k_pre = p_k
    lll1 = likelihood_SP(sigma2, sigma1_2, mu1, p_k, sigma2_2, index, wt, Y1 = Y1,
                         eta1_pre = eta1_pre,
                         eta2_pre = eta2_pre,
                         eta3_pre = eta3_pre,
                         k = k, X1 = X1, x = x, phi = 1/phi,
                         G = G, tttt = tttt)
    tmp = lll1$tmp
    G_ugk = foreach(gg = c(1:G), .combine = "rbind", .export = c("tmp", "k")) %do%{
      res = rep(0, k)
      if(max(tmp[gg,]) <= -700){
        res[which.max(tmp[gg,])] = 1
        res[which(tmp[gg, ] == -2000)] = 0
      }else{
        res = ((exp(tmp[gg,] ))/sum(exp(tmp[gg,])))
      }
      res
    }
    p_k = apply(G_ugk, 2, sum)/sum(apply(G_ugk, 2, sum))
    par = c(sigma1_2, mu1, sigma2, sigma2_2)
    result2 = nlminb(par, likelihood_f_sig1, lower = c(1e-10, -Inf, rep(1e-10, n_basis), rep(1e-10, n_basis)),
                     upper = rep(20, length(par)), control=list(iter.max = 1, trace=1, step.max = 1/(2^(j-1))))
    par = result2$par
    sigma1_2 = par[1]
    mu1 = par[2]
    sigma2 = par[(1:n_basis)+2]
    sigma2_2 = par[(1:n_basis) + 2 + n_basis]

    bb = proc.time()
    ll1 = likelihood_SP(sigma2, sigma1_2, mu1, p_k, sigma2_2, index, wt,
                        Y1 = Y1,
                        eta1_pre = eta1_pre,
                        eta2_pre = eta2_pre,
                        eta3_pre = eta3_pre,
                        k = k, X1 = X1, x = x, phi = 1/phi,
                        G = G, tttt = tttt)$ll
    proc.time() - bb
    index = list()
    wt = list()
    if(j == 1){
      n_step_pre <- nn
    }else{
      n_step_pre <- n_step
    }
    # n_step_pre <- n_step
    n_step <- round(n_step_pre/(2^(j-1)))
    if(n_step <= 100){
      n_step = 100
      j = j - 1
    }
    for (kk in c(1:k)){
      stat1 = apply(lll1$res_pre[[kk]], 2, sum)
      index[[kk]] = sample(x = c(1:(n_step_pre)), size = n_step, prob = (stat1 - min(stat1) + 1 )/sum((stat1 - min(stat1) + 1)))
      wt[[kk]] = (n_step * stat1/sum(stat1))[index[[kk]]]

    }
    ress = list(mu1 = mu1, sigma1_2 = sigma1_2, sigma2 = sigma2, sigma2_2 = sigma2_2, p_k = p_k)

  }
  result1 = list(mu1 = mu1, sigma1_2 = sigma1_2, sigma2 = sigma2, sigma2_2 = sigma2_2, p_k = p_k, phi = phi)

  return(list(data_use = data_use, result1 = result1))
}

