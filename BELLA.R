library(geepack)
library(coda)

######################### Functions ################################

## HMC function
HMC = function (log_U, grad_U, epsilon, L, q0){
  qq <- q0
  pp <- rnorm(length(qq), 0, 1)  
  current_p <- pp
  pp <- pp - epsilon * grad_U(qq) / 2
  
  for (i in 1:L){
    qq <- qq + epsilon * pp
    if (i != L) pp <- pp - epsilon * grad_U(qq)
  }
  pp <- pp - epsilon * grad_U(qq) / 2
  pp <- -pp
  
  current_U <- log_U(q0)
  current_K <- sum(current_p^2) / 2
  proposed_U <- log_U(qq)
  proposed_K <- sum(pp^2) / 2
  r <- exp(current_U - proposed_U + current_K - proposed_K)
  
  if (is.nan(r)) r <- 0
  
  if (runif(1) < r) {  
    th_new <- qq
    accepted <- TRUE
  } else {  
    th_new <- q0
    accepted <- FALSE
  }
  
  return(list(th = th_new, accepted = accepted))
}

### Functions about Ri matrix
Ri <- function(rho, ni){
  n  = length(ni)
  Ri = list()
  for(i in 1:n){
    Ri[[i]] = matrix(NA,ni[i],ni[i])
    Ri[[i]] = rho^(abs(row(Ri[[i]]) - col(Ri[[i]])))
  }
  return(Ri)
}

### Ri derivative
dRi <- function(rho, ni){
  n   = length(ni)
  dRi = list()
  for(i in 1:n){
    dRi[[i]] = matrix(NA,ni[i],ni[i])
    ind      = matrix(NA, nrow = ni[i], ncol = ni[i])
    indp     = abs(row(ind) - col(ind))
    dRi[[i]] = indp*rho^(indp-1)
  }
  return(dRi)
}

## Ri inverse
invRi <- function(rho,ni){
  n=length(ni)
  invRi=list()
  for(i in 1:n){
    invRi[[i]]=matrix(NA,ni[i],ni[i])
    if (ni[i] > 2) {
      diag(c(1, rep(1 + rho^2, ni[i] - 2), 1))
      invRi[[i]] = diag(c(1, rep(1+rho^2, ni[i]-2), 1))
      invRi[[i]][cbind(1:(ni[i]-1),2:ni[i])] = -rho
      invRi[[i]][cbind(2:ni[i],1:(ni[i]-1))] = -rho
      invRi[[i]] = (1/(1-as.numeric(rho)^2)) * invRi[[i]]
    } else {
      invRi[[i]] = matrix(c(1,-rho,-rho,1),2,2)
      invRi[[i]] = (1/(1-as.numeric(rho)^2)) * invRi[[i]]
    }
  } 
  return(invRi)
}

## Derivative of Ri inverse
dinvRi <- function(rho,ni){ 
  n=length(ni)
  dinvRi=list()
  for(i in 1:n){
    dinvRi[[i]]=matrix(NA,ni[i],ni[i])
    if (ni[i] > 2) {
      dinvRi[[i]] = diag(c(2*rho, rep(4*rho, ni[i]-2), 2*rho))
      dinvRi[[i]][cbind(1:(ni[i]-1),2:ni[i])] = -(1+rho^2)
      dinvRi[[i]][cbind(2:ni[i],1:(ni[i]-1))] = -(1+rho^2)
      dinvRi[[i]] = (1/((1-as.numeric(rho)^2)^2))*dinvRi[[i]]
    } else {
      dinvRi[[i]] = matrix(c(2*rho,-(1+rho^2),-(1+rho^2),4*rho),2,2)
      dinvRi[[i]] = (1/(1-as.numeric(rho)^2)^2) * dinvRi[[i]]
    } 
  }
  return(dinvRi)
}

## Second derivative of Ri inverse
ddinvRi <- function(rho, ni){
  n=length(ni)
  ddinvRi=list()
  for(i in 1:n){
    ddinvRi[[i]]=matrix(NA,ni[i],ni[i])
    if (ni[i] > 2) {
      ddinvRi[[i]] = diag(c((1+3*(rho^2)), rep((2)*(1+3*(rho^2)), ni[i]-2), (1+3*(rho^2))))
      ddinvRi[[i]][cbind(1:(ni[i]-1),2:ni[i])] = -rho*(3+rho^2)
      ddinvRi[[i]][cbind(2:ni[i],1:(ni[i]-1))] = -rho*(3+rho^2)
      ddinvRi[[i]] = (2/(1-as.numeric(rho)^2)^3) * ddinvRi[[i]]
    } else {
      ddinvRi[[i]] = matrix(c((1+3*(rho^2)),-rho*(3+rho^2),-rho*(3+rho^2),(1+3*(rho^2))),2,2)
      ddinvRi[[i]] = (2/(1-as.numeric(rho)^2)^3) * ddinvRi[[i]]
    } 
  }
  return(ddinvRi)
}

### Function for matrix G, with the estimating equations for beta, sigma2, and rho
mat_G=function(X, y, beta,sigma2,rho){
  
  if (!"id" %in% colnames(X)) stop("'X' must have 'id' column")
  
  ni      = as.vector(table(X$id)) 
  N       = sum(ni)
  n       = length(unique(X$id))     
  mat_G   = matrix(NA,(ncol(X)+1),n)
  Rii     = Ri(rho,ni)
  invRii  = invRi(rho,ni)
  dinvRii = dinvRi(rho,ni)
  i=1
  for(i in 1:n){
    iseq = X$id==i
    Xi   = X[iseq,-1]
    Xi   = apply(Xi, 2, as.numeric)
    
    #Estimating Equations
    eebeta = t(Xi)%*%invRii[[i]]%*%(y[iseq]-Xi%*%beta)
    eesig  = t((y[iseq]-Xi%*%beta))%*%invRii[[i]]%*%(y[iseq]-Xi%*%beta) - ni[i]*sigma2
    eerho  = (1/sigma2)*(t((y[iseq]-Xi%*%beta))%*%dinvRii[[i]]%*%(y[iseq]-Xi%*%beta)) - 2*rho*(ni[i]-1)/(1-rho^2)
    
    mat_G[,i] =  rbind(eebeta,eesig,eerho)
  }
  return(mat_G)
}

#### Function for Lagrange multiplier calculation
mylag <- function(G, lamb, gradtol=1e-7, decay_from=18){
  
  n = ncol(G)
  q = nrow(G)  
  
  if(!missing(lamb) ){
    lamb <- as.vector(lamb)
  }
  if(missing(lamb)) lamb <- rep(0,q)
  
  ## Pseudo log function L(lambda) - minimization over lambda
  llog <- function( z, eps ){
    ans <- z
    avoidNA <- !is.na(z)
    lo <- (z<eps) & avoidNA  
    ans[ lo  ] <- log(eps) - 1.5 + 2*z[lo]/eps - 0.5*(z[lo]/eps)^2
    ans[ !lo ] <- log( z[!lo] )
    return(ans)
  }  
  
  ## First derivative of L(lambda)
  llogp <- function( z, eps ){
    ans <- z
    avoidNA <- !is.na(z)   
    lo <- (z<eps) & avoidNA
    ans[ lo  ] <- 2.0/eps - z[lo]/eps^2
    ans[ !lo ] <- 1/z[!lo]
    return(ans)
  }
  
  ## Second derivative of L(lambda)    
  llogpp <- function( z, eps ){
    ans <- z
    avoidNA <- !is.na(z) 
    lo <- (z<eps) & avoidNA    
    ans[ lo  ] <- -1.0/eps^2
    ans[ !lo ] <- -1.0/z[!lo]^2
    return(ans)
  }
  
  
  ## Newton Raphson Iteractions 
  k = 0
  gam = 1  
  z = c(1 + t(lamb) %*% G)  
  dL1 = - apply(sweep(G, 2, llogp(z,1/n), "*"),1,sum) ##Calculations in BELLA document (p.17)
  d1size = sum(abs(dL1))
  
  while (TRUE){
    if (d1size < gradtol){
      break
    }
    dL2i=list()
    for(i in 1:n){
      dL2i[[i]]=-(G[,i]%*%t(G[,i]))*llogpp(z,1/n)[i]  ##Calculations in BELLA document (p.17)
    }
    dL2 = Reduce("+", dL2i)
    invdL2 = chol2inv(chol(dL2))   
    delta = gam * invdL2 %*% dL1
    lamb_cand = lamb - delta
    
    while (any(lamb_cand != lamb)){
      z_cand = c( 1 + t(lamb_cand) %*% G)                  
      dL1_cand = - apply(sweep(G, 2, llogp(z_cand,1/n), "*"),1,sum)
      d1size_cand = sum(abs(dL1_cand))
      
      if(d1size_cand < d1size){
        break
      }
      delta=delta/2
      lamb_cand= lamb - delta 
    }
    
    if(all(lamb_cand == lamb)){
      break
    }
    
    lamb = lamb_cand
    z = z_cand; dL1 = dL1_cand; d1size = d1size_cand
    
    k = k + 1
    if (k > decay_from){
      gam = 1/sqrt(k)
    }
  }  
  wts=1/(n*z)
  return(list(lamb=c(lamb), k=k, wts=wts))
}  


### Estimating beta, rho and sigma2
BELLA_thin <- function(X,y,nsim, nwarm, thin,
                  seed=sample.int(.Machine$integer.max, 1),
                  epsilon,L, a, b, gradtol=1e-7){
  time1 <- Sys.time()
  set.seed(seed)
  if (!"id" %in% colnames(X)) stop("'X' deve conter a coluna 'id'") 
  if (missing(epsilon)) epsilon = .8
  if (missing(thin))    thin    = 1
  
  X     <- as.data.frame(X)
  x     <- as.data.frame(X[,-1])
  y     <- as.vector(y)
  p     <- ncol(X)-1
  n     <- length(unique(X$id))      
  ni    <- as.vector(table(X$id))
  total <- floor((nsim - nwarm) / thin)
  theta.array       <- array(NA, dim = c(total,p+2))
  G.array           <- array(NA, dim = c(p+2,n,total))
  wts.array         <- array(NA, dim=c(total,n))
  acceptance_theta  <- rep(NA, total - 1)
  
  modelg = geeglm(y ~ ., id = X$id, data = data.frame(x[,-1]), family = gaussian, corstr = "ar1")
  theta.current <- c(modelg$coefficients, as.numeric(summary(modelg)$dispersion[1]), as.numeric(summary(modelg)$corr[1]))
  
  sample_id <- 1
  
  for(k in 1:nsim){
    log_U <- function(theta){
      matG = mat_G(X, y, theta[1:p],theta[p+1],theta[p+2])
      lam <-  mylag(matG)$lamb
      arg1 <- 0
      
      for(i in 1:n){
        arg1 <- arg1 + log(1 + t(lam) %*% matG[, i])
      }
      return(arg1 + ((0.5*t(theta[1:p])%*%solve(diag(p))%*%theta[1:p])/theta[p+1]) + b/theta[p+1] - (a+1)*log(theta[p+1]))
    }
    
    grad_U <- function(theta){
      matG = mat_G(X, y, theta[1:p],theta[p+1], theta[p+2])
      lam <-  mylag(matG)$lamb
      Rii=Ri(theta[p+2],ni)
      invRii=invRi(theta[p+2],ni)
      dinvRii=dinvRi(theta[p+2],ni)
      ddinvRii=ddinvRi(theta[p+2],ni)
      
      arg1 <- 0
      arg2 <- 0
      arg3 <- 0
      
      for(i in 1:n){
        iseq = X$id==i
        Xi   = X[iseq,-1]
        Xi   = apply(Xi, 2, as.numeric)
        
        ##Derivative of estimating equations 
        deeb_betai = (-1) * (t(Xi)%*%invRii[[i]]%*%Xi)
        deeb_sigi  = as.matrix(rep(0,p),col=1)
        deeb_rhoi  = as.numeric(t(Xi)%*%dinvRii[[i]]%*%(y[iseq]-Xi%*%theta[1:p]))
        dees_betai = 2 * (t(y[iseq]-Xi%*%theta[1:p])%*%invRii[[i]] %*% (-Xi))
        dees_sigi  = -ni[i]
        dees_rhoi  = as.numeric((t(y[iseq]-Xi%*%theta[1:p])%*%dinvRii[[i]]%*%(y[iseq]-Xi%*%theta[1:p])))
        deer_betai = (2/(theta[p+1]^2)) * (t((y[iseq]-Xi%*%theta[1:p])) %*% dinvRii[[i]] %*%(-Xi))
        deer_sigi  = as.numeric((-1/(theta[p+1])^2) * (t(y[iseq]-Xi%*%theta[1:p])%*%dinvRii[[i]]%*%(y[iseq]-Xi%*%theta[1:p]))) 
        deer_rhoi  = as.numeric((1/theta[p+1]) * t((y[iseq]-Xi%*%theta[1:p])) %*% ddinvRii[[i]] %*% ((y[iseq]-Xi%*%theta[1:p])) - 2*(ni[i]-1)*(1+theta[p+2]^2)/((1-theta[p+2]^2)^2))
        
        ##Summation
        arg1 <- arg1 + ((t(lam)%*%rbind(deeb_betai,dees_betai,deer_betai))/as.numeric(1 + t(lam) %*% matG[,i])) 
        arg2 <- arg2 + ((t(lam)%*%rbind(deeb_sigi,dees_sigi,deer_sigi))/as.numeric(1 + t(lam) %*% matG[,i])) 
        arg3 <- arg3 + ((t(lam)%*%c(deeb_rhoi, dees_rhoi, deer_rhoi))/as.numeric(1 + t(lam) %*% matG[,i])) 
      }
      
      arg1 <- arg1 + (1/theta[p+1])*t(theta[1:p])%*%solve(diag(p))
      arg2 <- arg2 - (1/(theta[p+1])^2)*(0.5*t(theta[1:p])%*%solve(diag(p))%*%theta[1:p] + b) - (a+1)/theta[p+1]
      
      return(c(arg1, arg2, arg3))
    }
    
    
    next.sample_theta = HMC(log_U = log_U, 
                            grad_U = grad_U, 
                            epsilon, 
                            L, q0 = theta.current)
    
    theta.current = next.sample_theta$th
    acceptance_theta[k-1] <- next.sample_theta$accepted
    
    ## Thinning
    if (k > nwarm && (k - nwarm) %% thin == 0) {
      theta.array[sample_id, ] <- theta.current
      G.array[,, sample_id] <- mat_G(X, y, theta.current[1:p], theta.current[p+1], theta.current[p+2])
      wts.array[sample_id, ] <- mylag(G = G.array[,, sample_id])$wts
      sample_id <- sample_id + 1
    }
    
    
    acceptance.rate_theta <- sum(acceptance_theta[-(1:nwarm)]) / (nsim - nwarm - 1)
    
    print(paste("I:", k, "Theta:", paste(format(theta.current, digits = 4), collapse = " ")))
  
  }
  time2 <- Sys.time()
  return(list(beta = theta.array[,1:p], sigma2 = theta.array[,(p+1)],
              rho = theta.array[,(p+2)], wts = wts.array,
              acceptance.theta=acceptance.rate_theta, 
              epsilon = epsilon,
              elapsedTime = difftime(time2, time1, units = "mins")))
}

BELLA_chain_thin <- function(X,y,nsim, nwarm, nchain, thin,
                        seed=sample.int(.Machine$integer.max,1),
                        epsilon, L, a, b, gradtol=1e-7) {
  time1 <- Sys.time()
  if (missing(thin))    thin    = 1
  X     <- as.data.frame(X)
  p     <- ncol(X) - 1
  ndraw <- floor((nsim - nwarm) / thin)
  result.chain = list()
  beta_list = lapply(1:nchain, function(x) matrix(NA, ndraw, p))
  sigma_mat = matrix(NA, ndraw, nchain) 
  rho_mat   = matrix(NA, ndraw, nchain) 
  wts_list   = lapply(1:nchain, function(x) matrix(NA, ndraw, n))
  acceptance =  NULL
  
  for (ii in 1:nchain){
    result.chain[[ii]]= BELLA_thin(X,y,nsim=nsim,nwarm=nwarm,seed=seed+ii,epsilon=epsilon,L=L,a=a,b=b, thin=thin)
    beta_list[[ii]] = result.chain[[ii]]$beta
    rho_mat[,ii]    = result.chain[[ii]]$rho
    sigma_mat[,ii]  = result.chain[[ii]]$sigma
    wts_list[[ii]]  = result.chain[[ii]]$wts
    acceptance = c(acceptance,result.chain[[ii]]$acceptance.theta)
  }
  time2 <- Sys.time()
  return(list(beta = beta_list, rho=rho_mat, 
              sigma = sigma_mat, wts = wts_list ,
              acceptance.rate=acceptance,
              epsilon = epsilon,
              elapsedTime = difftime(time2, time1, units = "mins")))
  
}





