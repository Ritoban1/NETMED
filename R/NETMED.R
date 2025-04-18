#' NETMED: Network-Based Mediation Analysis
#'
#' This function performs network-based mediation analysis using adjacency matrices,
#' mediators, exposures, and outcomes.
#'
#' @param A Numeric vector of exposures.
#' @param M Numeric vector of mediators.
#' @param Y Numeric vector of outcomes.
#' @param E Symmetric adjacency matrix (entries between 0 and 1; diagonal entries must be 1).
#' @param C Optional numeric matrix of covariates.
#' @param C_names Optional character vector of covariate names.
#' @return List containing the results of the network mediation analysis.
#' @import Matrix pracma MASS
#' @export

NETMED<-function(A,M,Y,E,C=NULL,C_names=NULL){
  require(Matrix)
  require(pracma)
  require(MASS)
  trace<-function(mat){
    return(sum(diag(mat)))
  }
  
  if(!is.vector(A)||!is.numeric(A))
    stop("A, vector of exposures should be a numeric vector")
  
  if(!is.vector(M)||!is.numeric(M))
    stop("M, vector of mediators should be a numeric vector")
  
  if(!is.vector(Y)||!is.numeric(Y))
    stop("Y, vector of outcomes should be a numeric vector")
  
  if(!is.matrix(E))
    stop("The adjacency matrix E should be a symmetric matrix where
         the entries are between 0 and 1 and the diagonal entries are 1.")
  
  p=0
  if (!is.null(C)) {
    if (!is.matrix(C)) {
      stop("C should be a matrix")
    }
    if(!is.vector(C_names)||!is.character(C_names))
      stop("C_names should be a character vector")
    p=ncol(C)
  }
  colnames(C)=C_names
  n=nrow(E)
  E_rowsum=apply(E,2,sum)
  
  ## tilde E
  E_t=E
  for(i in 1:n){
    E_t[i,]=E[i,]/(E_rowsum[i]-E[i,i])
    E_t[i,i]=0
  }
  mle_fun_m<-function(par){
    par[4+2*p]=exp(par[4+2*p])
    par[5+2*p]=exp(par[5+2*p])
    X=par[5+2*p]*E%*%E+par[4+2*p]*diag(1,nrow=nrow(E))
    E_inv=ginv(X)
    
    det=determinant(X,logarithm = TRUE)
    residuals=(cbind(1,A,E_t%*%A,C,
                     E_t%*%C)%*%par[c(1:(3+2*p))]-M)
    
    fun1=-0.5*det$modulus
    fun2=-0.5*t(residuals)%*%
      E_inv %*% (residuals)
    
    return(-(as.numeric(fun1+fun2)))
    
  }
  
  mle_fun_y<-function(par){
    par[6+2*p]=exp(par[6+2*p])
    par[7+2*p]=exp(par[7+2*p])
    X=par[7+2*p]*E%*%E+par[6+2*p]*diag(1,nrow=nrow(E))
    E_inv=ginv(X)
    
    det=determinant(X,logarithm = TRUE)
    residuals=(cbind(1,A,E_t%*%A,M,E_t%*%M,C,
                     E_t%*%C)%*%par[c(1:(5+2*p))]-Y)
    
    fun1=-0.5*det$modulus
    
    fun2=-0.5*t(residuals)%*%
      E_inv %*% (residuals)
    
    return(-(as.numeric(fun1+fun2)))
  }
  
  grad_fun_mle_m<-function(par){
    sigma_m=exp(par[4+2*p])
    sigma_b_m=exp(par[5+2*p])
    gamma=par[c(1:(3+2*p))]
    
    Sigma_m=sigma_m*diag(1,nrow=nrow(E))+sigma_b_m*E%*%E
    
    X_m=cbind(1,A,E_t%*%A,C,
              E_t%*%C)
    Sigma_m_inv=solve(Sigma_m)
    W_m=M-X_m%*%gamma
    grad_gamma=(t(X_m)%*%Sigma_m_inv%*%X_m)%*%gamma-(t(X_m)%*%Sigma_m_inv%*%M)
    grad_sigma_m= (0.5*trace(Sigma_m_inv) - 0.5*t(W_m)%*%(Sigma_m_inv %*% Sigma_m_inv)%*%W_m)*sigma_m
    grad_sigma_b_m= (0.5*trace(Sigma_m_inv %*% E %*% E) - 0.5*t(W_m)%*%(Sigma_m_inv
                                                                        %*% E %*% E %*%
                                                                          Sigma_m_inv)%*%W_m)*sigma_b_m
    c(grad_gamma,grad_sigma_m,grad_sigma_b_m)
    
  }
  grad_fun_mle_y<-function(par){
    sigma_y=exp(par[6+2*p])
    sigma_b_y=exp(par[7+2*p])
    beta=par[c(1:(5+2*p))]
    
    Sigma_y=sigma_y*diag(1,nrow=nrow(E))+sigma_b_y*E%*%E
    X_y=cbind(1,A,E_t%*%A,M,E_t%*%M,C,
              E_t%*%C)
    
    Sigma_y_inv=solve(Sigma_y)
    W_y=Y-X_y%*%beta
    
    grad_beta=(t(X_y)%*%Sigma_y_inv%*%X_y)%*%beta-
      (t(X_y)%*%Sigma_y_inv%*%Y)
    
    grad_sigma_y= (0.5*trace(Sigma_y_inv) - 0.5*t(W_y)%*%
                     (Sigma_y_inv %*% Sigma_y_inv)%*%W_y)*sigma_y
    
    grad_sigma_b_y= (0.5*trace(Sigma_y_inv %*% E %*% E) -
                       0.5*t(W_y)%*%
                       (Sigma_y_inv %*% E %*% E %*% Sigma_y_inv)%*%W_y)*sigma_b_y
    
    c(grad_beta,grad_sigma_y,grad_sigma_b_y)
    
  }
  
  variance_mle<-function(par_y,par_m){
    var<-matrix(0,nrow=(length(par_y)+length(par_m)),ncol=(length(par_y)+length(par_m)))
    
    ## y model
    sigma_y=par_y[6+2*p]
    sigma_b_y=par_y[7+2*p]
    beta_1=par_y[4]
    beta_2=par_y[5]
    beta_3=c(par_y[1],par_y[2],par_y[3],par_y[c(6:(5+2*p))])
    
    Sigma_y=sigma_y^2*diag(1,nrow=nrow(E))+sigma_b_y^2*E%*%E
    
    ## m model
    sigma_m=par_m[4+2*p]
    sigma_b_m=par_m[5+2*p]
    
    Sigma_m=sigma_m^2*diag(1,nrow=nrow(E))+sigma_b_m^2*E%*%E
    gamma=c(par_m[c(1:(3+2*p))])
    NR=cbind(1,A,E_t%*%A,C,
             E_t%*%C)
    
    Delta=matrix(0,nrow=2*length(A),ncol=2*length(A))
    Delta[c(1:length(A)),c(1:length(A))]=diag(1,nrow=length(A))
    Delta[(length(A)+1):(2*length(A)),(length(A)+1):(2*length(A))]=
      diag(1,nrow=length(A))
    Delta[c(1:length(A)),(length(A)+1):(2*length(A))]=
      -beta_1*diag(1,nrow=length(A))-
      beta_2*E_t
    
    D=c(Y,M)
    W=Delta%*%D-c(NR%*%beta_3,NR%*%gamma)
    
    mat_1<-function(A){
      mat=matrix(0,nrow=(2*nrow(A)),ncol=(2*nrow(A)))
      mat[1:nrow(A),1:nrow(A)]=A
      return(mat)
    }
    mat_2<-function(A){
      mat=matrix(0,nrow=(2*nrow(A)),ncol=(2*nrow(A)))
      mat[1:nrow(A),(nrow(A)+1):(2*nrow(A))]=A
      return(mat)
    }
    mat_4<-function(A){
      mat=matrix(0,nrow=(2*nrow(A)),ncol=(2*nrow(A)))
      mat[(nrow(A)+1):(2*nrow(A)),(nrow(A)+1):(2*nrow(A))]=A
      return(mat)
    }
    Sigma_y_inv=solve(Sigma_y)
    Sigma_m_inv=solve(Sigma_m)
    ## variance
    O<-matrix(0,nrow=(length(par_y)+length(par_m)),ncol=(length(par_y)+length(par_m)))
    
    O[1,1]=t(M)%*%Sigma_y_inv%*%M ##beta_1
    O[2,2]=t(E_t%*%M)%*%Sigma_y_inv%*%(E_t%*%M) ##beta_2
    O[1,2]=O[2,1]=t(M)%*%Sigma_y_inv%*%(E_t%*%M) ##beta_1 and beta_2
    
    O[c(3:(5+2*p)),(3:(5+2*p))]=t(NR)%*%Sigma_y_inv%*%NR ## beta_3
    O[1,c(3:(5+2*p))]=O[c(3:(5+2*p)),1]=as.numeric(t(M)%*%Sigma_y_inv%*%NR) ## beta_1 and beta_3
    O[2,c(3:(5+2*p))]=O[c(3:(5+2*p)),2]=as.numeric(t(E_t%*%M)%*%Sigma_y_inv%*%NR) ## beta_2 and beta_3
    
    O[(6+2*p),(6+2*p)]=-1/2 * trace(Sigma_y_inv%*%Sigma_y_inv) +
      t(W)%*%(mat_1(Sigma_y_inv%*%Sigma_y_inv%*%Sigma_y_inv))%*% W  ## sigma_y
    O[(7+2*p),(7+2*p)]=-1/2 * trace(Sigma_y_inv%*%E%*%E%*%Sigma_y_inv%*%E%*%E) +
      t(W)%*%(mat_1(Sigma_y_inv%*%E%*%E%*%Sigma_y_inv%*%E%*%E%*%Sigma_y_inv))%*% W ## sigma_b_y
    O[(6+2*p),(7+2*p)]=O[(7+2*p),(6+2*p)]=-1/2 * trace(Sigma_y_inv%*%E%*%E%*%Sigma_y_inv) +
      t(W)%*%(mat_1(Sigma_y_inv%*%Sigma_y_inv%*%E%*%E%*%Sigma_y_inv))%*%W ## sigma_y and sigma_b_y
    O[1,(6+2*p)]=O[(6+2*p),1]=t(Delta%*%D)%*%mat_2(Sigma_y_inv%*%Sigma_y_inv)%*%D-
      t(c(NR%*%beta_3,NR%*%gamma))%*%mat_2(Sigma_y_inv%*%Sigma_y_inv)%*%D ## beta_1 and sigma_y
    O[1,(7+2*p)]=O[(7+2*p),1]=t(Delta%*%D)%*%mat_2(Sigma_y_inv%*%E%*%E%*%Sigma_y_inv)%*%D-
      t(c(NR%*%beta_3,NR%*%gamma))%*%mat_2(Sigma_y_inv%*%E%*%E%*%Sigma_y_inv)%*%D ## beta_1 and sigma_b_y
    O[2,(6+2*p)]=O[(6+2*p),2]=t(Delta%*%D)%*%mat_2(Sigma_y_inv%*%Sigma_y_inv%*% E_t)%*%D-
      t(c(NR%*%beta_3,NR%*%gamma))%*%mat_2(Sigma_y_inv%*%Sigma_y_inv%*% E_t)%*%D ## beta_2 and sigma_y
    O[2,(7+2*p)]=O[(7+2*p),2]=t(Delta%*%D)%*%mat_2(Sigma_y_inv%*%E%*%E%*%Sigma_y_inv%*% E_t)%*%D-
      t(c(NR%*%beta_3,NR%*%gamma))%*%mat_2(Sigma_y_inv%*%E%*%E%*%Sigma_y_inv%*% E_t)%*%D
    O[c(3:(5+2*p)),(6+2*p)]=O[(6+2*p),c(3:(5+2*p))]=t(rbind(NR,matrix(0,nrow=nrow(NR),ncol=ncol(NR))))%*%
      mat_1(Sigma_y_inv%*%Sigma_y_inv)%*%c(NR%*%beta_3,NR%*%gamma)-
      t(rbind(NR,matrix(0,nrow=nrow(NR),ncol=ncol(NR))))%*%
      mat_1(Sigma_y_inv%*%Sigma_y_inv)%*%Delta%*%D ## beta_3 and sigma_y
    O[c(3:(5+2*p)),(7+2*p)]=O[(7+2*p),c(3:(5+2*p))]=t(rbind(NR,matrix(0,nrow=nrow(NR),ncol=ncol(NR))))%*%
      mat_1(Sigma_y_inv%*%E%*%E%*%Sigma_y_inv)%*%c(NR%*%beta_3,NR%*%gamma)-
      t(rbind(NR,matrix(0,nrow=nrow(NR),ncol=ncol(NR))))%*%
      mat_1(Sigma_y_inv%*%E%*%E%*%Sigma_y_inv)%*%Delta%*%D ## beta_3 and sigma_b_y
    
    O[c((8+2*p):(10+4*p)),c((8+2*p):(10+4*p))]=t(NR)%*%Sigma_m_inv%*%NR ## gamma
    O[(11+4*p),(11+4*p)]=-1/2 * trace(Sigma_m_inv%*%Sigma_m_inv) +
      t(W)%*%(mat_4(Sigma_m_inv%*%Sigma_m_inv%*%Sigma_m_inv))%*% W  ## sigma_m
    O[(12+4*p),(12+4*p)]=-1/2 * trace(Sigma_m_inv%*%E%*%E%*%Sigma_m_inv%*%E%*%E) +
      t(W)%*%(mat_4(Sigma_m_inv%*%E%*%E%*%Sigma_m_inv%*%E%*%E%*%Sigma_m_inv))%*% W ## sigma_b_m
    O[(11+4*p),(12+4*p)]=O[(12+4*p),(11+4*p)]=
      -1/2 * trace(Sigma_m_inv%*%E%*%E%*%Sigma_m_inv) +
      t(W)%*%(mat_4(Sigma_m_inv%*%Sigma_m_inv%*%E%*%E%*%Sigma_m_inv))%*%W ## sigma_m and sigma_b_m
    
    O[c((8+2*p):(10+4*p)),(11+4*p)]= O[(11+4*p),c((8+2*p):(10+4*p))]=
      t(rbind(matrix(0,nrow=nrow(NR),ncol=ncol(NR)),NR))%*%
      mat_4(Sigma_m_inv%*%Sigma_m_inv)%*%c(NR%*%beta_3,NR%*%gamma)-
      t(rbind(matrix(0,nrow=nrow(NR),ncol=ncol(NR)),NR))%*%
      mat_4(Sigma_m_inv%*%Sigma_m_inv)%*%Delta%*%D ## gamma and sigma_y
    O[c((8+2*p):(10+4*p)),(12+4*p)]= O[(12+4*p),c((8+2*p):(10+4*p))]=
      t(rbind(matrix(0,nrow=nrow(NR),ncol=ncol(NR)),NR))%*%
      mat_4(Sigma_m_inv%*%E%*%E%*%Sigma_m_inv)%*%c(NR%*%beta_3,NR%*%gamma)-
      t(rbind(matrix(0,nrow=nrow(NR),ncol=ncol(NR)),NR))%*%
      mat_4(Sigma_m_inv%*%E%*%E%*%Sigma_m_inv)%*%Delta%*%D ## gamma and sigma_b_y
    
    x=rep(0,times=(13+4*p))
    inv_O=solve(O)
    x[c((8+2*p):(12+4*p))]=diag(inv_O)[c((8+2*p):(12+4*p))]
    x[1]=diag(inv_O)[3]
    x[2]=diag(inv_O)[4]
    x[3]=diag(inv_O)[5]
    x[4]=diag(inv_O)[1]
    x[5]=diag(inv_O)[2]
    x[c(6:(7+2*p))]=diag(inv_O)[c(6:(7+2*p))]
    
    x[(13+4*p)]=inv_O[c(9+2*p),(10+2*p)]
    return(x)
  }
  
  interaction_terms <- sapply(C_names, function(C) paste("E_t %*%", C))
  
  formula_parts <- c("Y ~ A", "E_t %*% A", "M", "E_t %*% M",
                     paste(C_names, collapse = " + "), paste(interaction_terms, collapse = " + "))
  formula_str <- paste(formula_parts, collapse = " + ")
  
  model_formula <- as.formula(formula_str)
  Y_m <- lm(model_formula)
  sum_y=summary(Y_m)
  # Exclude M and its interaction from the model parts
  formula_parts_M <- c("M ~ A", "I(E_t %*% A)",
                       paste(C_names, collapse = " + "), paste(interaction_terms, collapse = " + "))
  formula_str_M <- paste(formula_parts_M, collapse = " + ")
  model_formula_M <- as.formula(formula_str_M)
  M_m <- lm(model_formula_M)
  sum_m=summary(M_m)
  try({
    b_y_sig=runif(1,0,log(var(Y)))
    model_y=optim(par=c(as.numeric(coef(Y_m)),
                        log(var(Y)),b_y_sig),mle_fun_y,
                  gr= grad_fun_mle_y,method="BFGS",
                  control=list(trace=TRUE))
  })
  
  betahat_y=model_y$par[c(1:(5+2*p))]
  sigma_y_hat=sqrt(exp(model_y$par[(6+2*p)]))
  sigma_b_y_hat=sqrt(exp(model_y$par[(7+2*p)]))
  
  mle_y=c(betahat_y,sigma_y_hat,sigma_b_y_hat)
  
  ## MLE of M model
  try({
    b_m_sig=runif(1,0,log(var(M)))
    model_m=optim(par=c(as.numeric(coef(M_m)),log(var(M)),b_m_sig),
                  mle_fun_m,gr=
                    grad_fun_mle_m,method="BFGS",
                  control=list(trace=TRUE))
    
  })
  gammahat_m=model_m$par[c(1:(3+2*p))]
  sigma_m_hat=sqrt(exp(model_m$par[(4+2*p)]))
  sigma_b_m_hat=sqrt(exp(model_m$par[(5+2*p)]))
  mle_m=c(gammahat_m,sigma_m_hat,sigma_b_m_hat)
  mle=c(mle_y,mle_m)
  var_par=variance_mle(mle_y,mle_m)
  
  var_par[which(var_par[c(1:(5+2*p))]<0)]=
    sum_y$coefficients[,2][which(var_par[c(1:(5+2*p))]<0)]
  
  var_par[which(var_par[c((8+2*p):(10+4*p))]<0)]=
    sum_m$coefficients[,2][(which(var_par[c((8+2*p):(10+4*p))]<0)-(7+2*p))]
  
  ####
  ## E1
  E1_est=mle[2]
  E1_var=var_par[2]
  
  #E2
  E2_est=mle[4]*mle[(9+2*p)]
  E2_var=var_par[4]*mle[(9+2*p)]^2+mle[4]^2*var_par[(9+2*p)]
  
  
  #E3
  ## E3
  
  N <- rowSums(E) - 1
  
  
  sum_3 <- 0
  
  
  for (i in 1:nrow(E)) {
    sum_3 <- sum_3 + sum(E[i, -i] / N[-i]) / N[i]
  }
  
  
  E3_est=mle[5]*mle[(10+2*p)]*1/nrow(E) * sum_3
  
  E3_var=(var_par[5]*mle[(10+2*p)]^2+mle[5]^2*var_par[(10+2*p)])*(1/nrow(E) * sum_3)^2
  
  ## E4
  E4_est=mle[3]
  
  E4_var=var_par[3]
  
  #E5
  E5_est=mle[4]*mle[(12+2*p)]
  
  
  E5_var=var_par[4]*mle[(10+2*p)]^2+mle[4]^2*var_par[(10+2*p)]
  
  ### E6
  inv_N <- 1 / N
  
  sum_6 <- 0
  
  
  for (i in 1:nrow(E)) {
    for (j in 1:nrow(E)) {
      if (j != i) {
        valid_k <- (1:nrow(E) != i) & (1:nrow(E) != j) & (E[i, ] == 1)
        if (any(valid_k)) {
          sum_6 <- sum_6 + sum(E[i, j] * E[j, valid_k] * inv_N[j]) * inv_N[i]
        }
      }
    }
  }
  E6_est=mle[5]*mle[(11+2*p)]+mle[5]*mle[(12+2*p)]*sum_6*1/nrow(E)
  
  ## var
  E6_var=(mle[(9+2*p)]+mle[(10+2*p)]*sum_6*1/nrow(E))^2*
    var_par[5]+mle[5]^2*(var_par[(9+2*p)]+var_par[(10+2*p)]*(sum_6*1/nrow(E))^2) +
    2*mle[5]^2*sum_6*1/nrow(E)*var_par[(13+4*p)]
  
  res <- data.frame(
    Est = c(E1_est, E2_est, E3_est, E4_est, E5_est, E6_est),
    Var = c(E1_var, E2_var, E3_var, E4_var, E5_var, E6_var)
  )
  
  
  return(list(mle_est=mle,var_mle=var_par,Effects_res=res))
  
}



