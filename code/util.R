
update_q <- function(x, s, mu, sigma2, init=list(NULL), control=list(maxiter=100, delta=1, tol=1e-8, lwr=-6, upr=10)){
  n <- length(x)
  K <- length(sigma2)
  maxiter <- control$maxiter
  delta <- control$delta
  tol <- control$tol
  lwr <- control$lwr
  upr <- control$upr
  m <- init$m
  v2 <- init$v2
  
  if(is.null(m)){
    m <- matrix(0, nrow=n, ncol=K)
  }
  
  if(is.null(v2)){
    v2 <- matrix(1, nrow=n, ncol=K)
  }
  
  sigma2_mat <- rep(1, n) %*% t(sigma2)
  x_mat <- x %*% t(rep(1, K))
  s_mat <- s %*% t(rep(1, K))
  a <-  diag(s)%*%exp(mu + m + v2/2)
  
  for(iter in 1:maxiter){
    v2_new <- 1/(1/sigma2_mat + a)
    
    m_new <- m - delta*v2_new*(a - x_mat + m/sigma2_mat)
    
    idx.lwr <- (m_new + log(s_mat) < lwr - mu)
    if(sum(idx.lwr) !=0){
      m_new[idx.lwr] <- lwr - mu - log(s_mat[idx.lwr]) 
    }
    
    idx.upr <- (m_new + log(s_mat) > upr - mu)
    if(sum(idx.upr) !=0){
      m_new[idx.upr] <- upr - mu - log(s_mat[idx.upr]) 
    }    
    
    if(max(abs(m_new-m)) < tol & max(abs(v2_new - v2)) < tol) break
    
    m <- m_new
    v2 <- v2_new
    a <- diag(s)%*%exp(mu + m + v2/2)
  }
  
  ELBO <- x_mat*m - a - 0.5*log(sigma2_mat) + 0.5*log(v2) - 0.5*(m^2+v2)/sigma2_mat + 0.5
  
  return(list(mu=mu, sigma2=sigma2, m=m, v2=v2, a=a, ELBO=ELBO))
}


pois_ash <- function(x, s, sigma2, init=NULL, maxiter=100, tol=1e-6, verbose=FALSE){
  n <- length(x)
  K <- length(sigma2)
  
  mu <- init$mu
  if(is.null(mu)){
    mu <- log(sum(x)) - log(sum(s))
  }
  
  pi <- init$pi
  if(is.null(pi)){
    pi <- rep(1/K, K)
  }

  const <- sum(x*log(s)) - sum(lgamma(x+1))
  ELBOs <- c()
  
  for(iter in 1:maxiter){
    # update posterior mean m_ik and variance v2_ik
    res.q <- update_q(x, s, mu, sigma2)
    tmp.q <- diag(s)%*%exp(res.q$m + res.q$v2/2)
    
    # update posterior mean of z_ik
    ELBO.local <- res.q$ELBO
    ELBO.cen <- ELBO.local - apply(ELBO.local, 1, max)
    Ez <- exp(ELBO.cen)%*%diag(pi)
    Ez <- diag(1/rowSums(Ez)) %*% Ez
    Ez <- pmax(Ez, 1e-15)
    
    # compute overall ELBO
    pi_mat <- rep(1, n) %*% t(pi)
    ELBO.overall <- mu*sum(x) + sum(Ez*(log(pi_mat) + ELBO.local - log(Ez))) + const
    ELBOs <- c(ELBOs, ELBO.overall) 
    
    if(verbose){
      print("iter         ELBO")
      print(sprintf("%d:    %f", iter, ELBO.overall))
    }
    
    # update mu
    mu.new <- log(sum(x)) - log(sum(Ez*tmp.q)) 
    diff.mu <- mu.new - mu
    mu <- mu.new
    
    # update pi
    pi.new <- colMeans(Ez)
    pi.new <- pmax(pi.new, 1e-8)
    diff.pi <- pi.new - pi
    pi <- pi.new
    
    # if(abs(diff.mu) < tol & max(abs(diff.pi)) < tol) break
  }
  
  return(list(mu=mu, diff.mu=diff.mu, sigma2=sigma2, pi=pi, diff.pi=diff.pi, Ez=Ez, m=res.q$m, v2=res.q$v2, ELBO=ELBOs))
}


