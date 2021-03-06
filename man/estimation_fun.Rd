% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation_fun.R
\name{estimation_fun}
\alias{estimation_fun}
\title{Estimation for MAP data set:}
\usage{
estimation_fun(n_control = 10, n_treat = 10, n_rep = 3, x, Y1,
  k = 4, nn = 300, phi = NULL, type = NULL, tttt)
}
\arguments{
\item{n_control}{Number of time points in control group}

\item{n_treat}{Number of time points in treatment group}

\item{n_rep}{Number of replicates in both group}

\item{x}{Time structured design for the simulated data}

\item{Y1}{RNA-seq data set: G by (n_control + n_treat) * n_rep matrix}

\item{k}{Always 4}

\item{nn}{Number of QMC nodes used in likelihood estimation}

\item{phi}{the dispersion parameter}

\item{type}{dispersion parameter estimation type}

\item{tttt}{Time points used in the design}
}
\value{
A list of parameters that we need to use for DE analysis.\cr

data_use List includes the data information
\itemize{
\item Y1 RNA-seq data set: G by (n_control + n_treat) x n_rep matrix
\item start Starting point of the EM algorithm
\item k always 4
\item n_basis Number of basis function have been used to construct the design matrix
\item X1 Vector indicate control data (0) or treatment data (1)
\item x Design matrix been used for estimation
\item tttt Time points used in the design \cr
}

result1 List includes the Parameters estimations
\itemize{
\item mu1 Mean parameter for \eqn{\eta_{g1}}
\item sigma1_2 Variance parameter for \eqn{\eta_{g1}}
\item sigma2 Variance parameter for \eqn{\tau_{g}}
\item sigma2_2 Variance parameter for \eqn{\eta_{g2}}
\item p_k Proportion for each mixture component
\item aa Dispertion parameter
}
}
\description{
Estimation for MAP data set:
}
\details{
The vector of read counts for gene g, treatment group i, replicate j,
at time point \eqn{t,Y_{gij}(t)}, follows a Negative Binomial distribution
parameterized mean \eqn{\lambda_{gi}} and \eqn{\phi_g}, where
\eqn{E[Y_{gij}(t)] = \lambda_{gi}(t)}.
\eqn{\lambda_{gi}(t)} is further modeled as
\eqn{\lambda_{gi}(t) = S_{ij} \exp[\eta_{g1}I_{i = 2} + B'(t)\eta_{g2}I_{i = 2} + B'(t)\tau_{g}].}
We have \eqn{B'(t)} are design matrix, which is constructed by 2 orthorgonal polynomial bases.
\itemize{
\item t = 1,..., n_treat (or n_control if control group);
\item j = 1,..., n_rep;
\item g = 1,...,G; and
\item \eqn{[\eta_{g1}, \eta_{g2}, \tau_{g}]} ~ 4-component gausssian mixture model.
We used latented negative binomial model with EM algorithm to estimate the paramters of mixture model.
}
}
\examples{
library(matlib)
n_basis = 2
n_control = 10
n_treat   = 10
n_rep = 3
tt_treat  = c(1:n_treat)/n_treat
nt = length(tt_treat)
ind_t = sort(sample(c(1:nt), n_control))
tt = tt_treat[ind_t]
tttt = c(rep(tt, n_rep), rep(tt_treat, n_rep))
z = x = matrix(0, length(tttt), n_basis)
z[,1] = 1.224745*tttt
z[,2] = -0.7905694 + 2.371708*tttt^2
x[,1] = z[,1] - Proj(z[,1], rep(1, length(tttt)))
x[,2] = z[,2] - Proj(z[,2], rep(1, length(tttt))) - Proj(z[,2], x[,1])
Y1 = data_generation(G = 100,
                     n_control = n_control,
                     n_treat   = n_treat,
                     n_rep     = n_rep,
                     k_real = 4,
                     sigma2_r = rep(1, 2),
                     sigma1_2_r = 1,
                     sigma2_2_r = c(3,2),
                     mu1_r = 4,
                     phi_g_r = rep(1,100),
                     p_k_real = c(0.7, 0.1, 0.1, 0.1),
                     x = x)
aaa <- proc.time()
est_result <- estimation_fun(n_control = n_control,
                             n_treat = n_treat,
                             n_rep = n_rep,
                             x = x,
                             Y1 = Y1,
                             nn = 300,
                             k = 4,
                             phi = NULL,
                             type = 2,
                             tttt = tttt)
aaa1 <- proc.time()
aaa1 - aaa
#------------gaussian basis construction------------------
# method <- "gaussian"
# n_basis <- 2
# a = 1/4
# b = 2
# c = sqrt(a^2 + 2 * a * b)
# A = a + b + c
# B = b/A
# phi_cal = function(k, x){
#   Lambda_k =sqrt(2 * a / A) * B^k
#   1/(sqrt(sqrt(a/c) * 2^k * gamma(k+1))) *
#      exp(-(c-a) * x^2) * hermite(sqrt(2 * c) * x, k, prob = F)
#
# }
# z = do.call(cbind, lapply(c(1:n_basis), phi_cal, x = (tttt - mean(tttt))))
# x = matrix(0, length(tttt), n_basis)
# if(n_basis == 2){
#   x[,1] = z[,1] - Proj(z[,1], rep(1, length(tttt)))
#   x[,2] = z[,2] - Proj(z[,2], rep(1, length(tttt))) - Proj(z[,2], x[,1])
#
# }else{
#   x[,1] = z[,1] - Proj(z[,1], rep(1, length(tttt)))
#   x[,2] = z[,2] - Proj(z[,2], rep(1, length(tttt))) - Proj(z[,2], x[,1])
#   x[,3] = z[,3] - Proj(z[,3], rep(1, length(tttt))) - Proj(z[,3], x[,1]) - Proj(z[,3], x[,2])
# }
}
