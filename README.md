# MAPTest
## Prerequisites
You also need an R (>= 3.4.0) installation on your computer; run:

```
install.packages(c("foreach","MASS","mvtnorm","randtoolbox","stats","mclust","EQL","matlib","parallel"), dependencies=TRUE)
```
which should download, build and install all R tools needed.

## Build, install and test

Run R shell and load devtools library:
```
require(devtools)

install_github("meca7653/MAPTest")

library(MAPTest)
```

## Examples:

### Simulate Data

The Example_Data available on Biometrics website on Wiley Online Library can be generated from the following code.

```
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

p_k_real = c(0.5, 0.2, 0.2, 0.1)

set.seed(2019)
Y1 = data_generation(G = 1000,
                     n_control = n_control,
                     n_treat   = n_treat,
                     n_rep     = n_rep,
                     k_real = 4,
                     sigma2_r = rep(1, 2),
                     sigma1_2_r = 1,
                     sigma2_2_r = c(3,2),
                     mu1_r = 4,
                     phi_g_r = rep(1, 1000),
                     p_k_real = p_k_real,
                     x = x)
colnames(Y1) <- c(paste(paste("C", rep(c(1:n_control), n_rep), sep = "_"), rep(c(1:n_rep), each = n_control), sep = "_"),
                  paste(paste("T", rep(c(1:n_treat), n_rep), sep = "_"), rep(c(1:n_rep), each = n_treat), sep = "_"))
rownames(Y1) <- paste("Gene", c(1:1000))
save(Y1, file = "Example_Data.rda")
write.csv(Y1, file = "Example_Data.csv")
```

### Estimate

```
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
```
### MAPTest

```
G <- 1000

k_real <- 4

p_k_real <- c(0.5, 0.2, 0.2, 0.1)

dd = rep(c(0:(k_real-1)), p_k_real * G)

result = MAP_test(est_result = est_result, Type = c(1:6), dd = dd, nn = 300)

Summary_MAP(result)
```
