##############################################################
##                Functions used in this study              ##
##                          04/05/2018                      ##
##    1. create moment/gram tensors from data               ##
##    2. create cumulant tensors from data                  ##
##    3. HOEVD: hosvd for symmetric tensors                 ##
##    4. hosvd for gram matrix with large n                 ##
##       with reduction from n to p                         ##
##       useful to get U from USV'                          ##
##    5. robust tensor power method with                    ##
##       vectorization                                      ## 
##############################################################
library(rTensor)
#1. moment/gram tensors

# 1. create momentGram tensors --------------------------------------------
MGT3<- function (x) {
  n <- nrow (x)
  m <- ncol (x)
  mm <- 1 : m
  nn <- 1 : n
  c3 <- array (0, c (m, m, m))
  for (i in nn) {
    c3 <- c3 + outer (outer (x [i, ], x [i, ]), x [i,])
  }
  return (c3)
}


MGT4 <- function (x) {
  n <- nrow (x)
  m <- ncol (x)
  mm <- 1 : m
  nn <- 1 : n
  c4 <- array (0, c (m, m, m, m))
  for (i in nn) {
    c4 <- c4 + outer (outer (x [i, ], x [i, ]), outer (x [i, ], x [i, ]))
  }
  return (c4)
}


# 2. Create cumulan tensors -----------------------------------------------
third_cumulant<- function (x) {
  n <- nrow (x)
  m <- ncol (x)
  mm <- 1 : m
  nn <- 1 : n
  c3 <- array (0, c (m, m, m))
  for (i in nn) {
    c3 <- c3 + outer (outer (x [i, ], x [i, ]), x [i,])
  }
  return (c3)
}


four_cumulants <- function (x) {
  n <- nrow (x)
  m <- ncol (x)
  mm <- 1 : m
  nn <- 1 : n
  r2 <- crossprod (x) / n
  r4 <- array (0, c (m, m, m, m))
  for (i in nn) {
    r4 <- r4 + outer (outer (x [i, ], x [i, ]), outer (x [i, ], x [i, ]))
  }
  r4<-r4 / n
  
  c2<-r2 
  c4<-r4
  for (i in mm) for (j in mm) for (k in mm) for (l in mm)
  {
    s4 <- r4 [i, j, k, l]
    s22 <- r2 [i, j] * r2 [k, l] + r2 [i, k] * r2 [j, l] + r2 [j, k] * r2 [i, l]
    c4 [i, j, k, l] <- s4 -  s22 
  }
  return (c4)
}


# 3. HOEVD: hosvd for symmetric tensors -----------------------------------
hoevd = function(x,rank){
  tnsr = as.tensor(x)
  unfold = k_unfold(tnsr,m = 1)@data
  u = svd(unfold, nu = rank)$u
  ulist = list(u)
  ulist = rep(ulist,tnsr@num_modes)
  z = ttl(tnsr, lapply(ulist, t), ms = 1:tnsr@num_modes)@data
  return(list( u = u, z = z))
}



# 4. hosvd for gram tensors with large n ----------------------------------
Gram3_hosvd = function(y){
  n = dim(y)[1]
  p = dim(y)[2]
  
  svd_y = svd(y)
  u2 = svd_y$u
  v2 = svd_y$v
  s2 = diag(svd_y$d)
  
  vvv = MGT3(v2)
  middle = s2 %*% k_unfold(as.tensor(vvv),m=1)@data %*% (s2 %x% s2)
  middle_u = svd(middle)$u
  
  
  u3 = u2 %*% middle_u
  return(list (u = u3, middle = middle))
}

Gram4_hosvd = function(y){
  n = dim(y)[1]
  p = dim(y)[2]
  
  svd_y = svd(y)
  u2 = svd_y$u
  v2 = svd_y$v
  s2 = diag(svd_y$d)
  
  vvvv = MGT4(v2)
  middle = s2 %*% k_unfold(as.tensor(vvvv),m=1)@data %*% (s2 %x% s2 %x% s2)
  middle_u = svd(middle)$u
  
  u4 = u2 %*% middle_u
  return(list (u = u4, middle = middle))
}
# 5. Robust tensor power method (vectorization) ---------------------------
# power iteration for the first eigen vecotr
power_itr1 = function(theta,Y,N,order){
  n = nrow(Y)
  d = 1/n
  
  for(k in 1:N){
    w = (Y %*% theta)^(order-1)
    next_itr = as.vector(t(w) %*% Y)
    next_itr = d * next_itr
    
    thetai = next_itr/sqrt(sum(next_itr^2))
    err = sum((thetai - theta)^2)
    theta = thetai
    if( err <= 1e-6){break}
  }
  return(theta)
}

# power iteration for all following eigen vectors
power_itrk = function(theta,eigenv, eigenl,Y,N,order){
  n = nrow(Y)
  d = 1/n
  
  for(k in 1:N){
    w1 = (Y %*% theta)^(order-1)
    next_itr1 = as.vector(t(w1) %*% Y)
    
    w2 = (t(eigenv) %*% theta)^(order - 1)
    if(length(eigenl) == 1){next_itr2 = as.vector(eigenv %*% eigenl %*% w2)}
    if(length(eigenl) >1){next_itr2 = eigenv %*% diag(eigenl) %*% w2}
    
    next_itr = d * next_itr1 - next_itr2
    thetai = next_itr/sqrt(sum(next_itr^2))
    err = sum((thetai - theta)^2)
    theta = thetai
    if( err <= 1e-6){break}
  }
  return(theta)
}

# get first lambda
Est_PI1 = function(Y, L=10, N=10,order){
  p = ncol(Y)
  theta_list = list()
  for ( t in 1:L){
    vc = rnorm(p)
    theta_0 = vc/sqrt(sum(vc^2))
    theta_list[[t]] = power_itr1(theta = theta_0, Y = Y, N=N, order = order)
  }
  
  n = nrow(Y)
  d = 1/n
  lambda_list = NULL
  
  for(t in 1:L){
    a = (Y %*% theta_list[[t]])^order
    lambda_t = d * sum(a)
    lambda_list = c(lambda_list,lambda_t)
  }
  
  ind = which.max(abs(lambda_list))
  theta_tau = theta_list[[ind]]
  
  theta_hat = power_itr1(theta = theta_tau,Y = Y, N=N, order = order)
  
  a = (Y %*% theta_hat)^order
  lambda_hat = d * sum(a)
  
  result = list("theta_hat" = theta_hat, "lambda_hat" = lambda_hat)
  return(result)
} 

# get all following lambda
Est_PIk = function(Y,eigenv, eigenl,L=10, N=10,order){
  p = ncol(Y)
  theta_list = list()
  for ( t in 1:L){
    vc = rnorm(p)
    theta_0 = vc/sqrt(sum(vc^2))
    theta_list[[t]] = power_itrk(theta = theta_0, eigenv, eigenl, Y = Y, N=N, order = order)
  }
  
  n = nrow(Y)
  d = 1/n
  lambda_list = NULL
  
  for(t in 1:L){
    a =  sum((Y %*% theta_list[[t]])^order)
    if(length(eigenl) == 1){ b =  sum(eigenl * ((t(eigenv) %*% theta_list[[t]])^order))}
    if(length(eigenl) >1){b =  sum(diag(eigenl) %*% ((t(eigenv) %*% theta_list[[t]])^order))}
    lambda_t = d * a - b
    lambda_list = c(lambda_list,lambda_t)
  }
  
  ind = which.max(abs(lambda_list))
  theta_tau = theta_list[[ind]]
  theta_hat = power_itrk(theta = theta_tau,eigenv, eigenl,Y = Y, N=N, order = order)
  
  a =  sum((Y %*% theta_hat)^order)
  if(length(eigenl) == 1){ b =  sum(eigenl * ((t(eigenv) %*% theta_hat)^order))}
  if(length(eigenl) >1){b =  sum(diag(eigenl) %*% ((t(eigenv) %*% theta_hat)^order))}
  
  lambda_hat = d * a - b
  
  result = list("theta_hat" = theta_hat, "lambda_hat" = lambda_hat)
  return(result)
} 

Rob_TPM = function(Y,L = 10, N = 10,order, p = ncol(Y)){
  
  result.1 = Est_PI1(Y = Y, L = L, N = N, order = order)
  eigenv = matrix(result.1$theta_hat,nrow = length(result.1$theta_hat))
  eigenl = result.1$lambda_hat
  
  for( m in 2:p){
    result.m = Est_PIk(Y = Y,eigenv = eigenv, eigenl = eigenl, L = L, N = N, order = order)
    eigenv = cbind(eigenv, result.m$theta_hat)
    eigenl = c(eigenl, result.m$lambda_hat)
    print(m)
  }
  
  eigenv = eigenv[,order(abs(eigenl),decreasing = T)]
  eigenl = eigenl[order(abs(eigenl),decreasing = T)]
  
  
  result = list("eigenv" = eigenv, "eigenl" = eigenl)
}