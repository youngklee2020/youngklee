# R functions for LL method for density response 

# Function for obtaining candidate bandwidths
# X: observed X (matrix of size (sample size,number of covariates))
# cv: number of folds for cross-validation
# h_add and h_length: ncol(X)-dimensional vectors
# seq(min_h[j],min_h[j]+h_add[j],length=h_length[j]) will be the candidate set of j-th bandwidth
# for some small bandwidth min_h[j]
# resulting h_candidate_matrix is the matrix whose j-th column is the candidate set of j-th bandwidth
get_h_candidate=function(X,cv,h_add,h_length)
{
  n=nrow(X)
  d=ncol(X)
  s=sample(1:n)
  X=X[s,]
  folds=new_cut(n,cv)
  distance=matrix(,nrow=cv,ncol=d)
  for(k in 1:cv)
  {
    X_training=X[-which(folds==k),]
    for(j in 1:d)
    {
      X_j=sort(X_training[,j])
      distance[k,j]=max((X_j[3:length(X_j)]-X_j[1:(length(X_j)-2)])/2,X_j[2],1-X_j[length(X_j)-1])
    }
  }
  min_h=c()
  h_candidate_matrix=matrix(,nrow=max(h_length),ncol=d)
  for(j in 1:d)
  {
    min_h[j]=max(distance[,j])+0.001
    h_candidate_matrix[1:h_length[j],j]=seq(min_h[j],min_h[j]+h_add[j],length=h_length[j])
  }
  return(list(h_candidate_matrix=h_candidate_matrix,s=s))
}

# Function giving index of splited data for cross-validation
# n: sample size
# cv: number of folds for cross-validation
# only the last fold has a different number of observations
new_cut=function(n,cv)
{
  folds=c()
  size=floor(n/cv)
  size_last=n-(cv-1)*size
  for(k in 1:cv)
  {
    if(k<cv) folds[((k-1)*size+1):(k*size)]=k
    if(k==cv) folds[((cv-1)*size+1):n]=k
  }
  return(folds)
}

# Function for obtaining additional grid points for numerical integration
# x: X_test matrix taking values in [0,1]^d
# we use sort(X_test[,j],seq(0,1,length=add)) as the grid points for numerical integration for j-th integral
integration.grid=function(x,add)
{
  add_vector=seq(0,1,length=add)
  added_x=matrix(,add+nrow(x),ncol(x))
  for(j in 1:ncol(x))
  {
    added_x[,j]=sort(c(add_vector,x[,j]))
  }
  added_x
}

# Function for CBS
# X_training and Y_training are usual
# add: number of additional grid points for numerical integration
# Y_grid: common observed times for Y
# max_iteration: maximum iteration number for Bochner smooth bacfkfitting algorithm
# epsilon: convergence criterion for Bochner smooth backfitting algorithm
# max_cbs_iteration: maximum iteration number for CBS algorithm
# cv: number of folds for cross-validation
# h_add and h_length: ncol(X_training)-dimensional vectors
# seq(min_h[j],min_h[j]+h_add[j],length=h_length[j]) will be the candidate set of j-th bandwidth
# for some small bandwidth min_h[j]
CBS_density=function(X_training,Y_training,add=11,Y_grid=NULL,max_iteration=50,epsilon=10^-4
                     ,max_cbs_iteration=20,cv=10,h_add=NULL,h_length=NULL)
{
  n=nrow(X_training)
  d=ncol(X_training)
  T=ncol(Y_training)
  
  if(is.null(h_add)) h_add=rep(0.5,d)
  if(is.null(h_length)) h_length=as.integer(rep(101,d))
  else h_length=as.integer(h_length)
  if(is.null(Y_grid)) Y_grid=as.double(seq(0,1,length=T))
  else Y_grid=as.double(Y_grid)
  add=as.integer(add)
  max_iteration=as.integer(max_iteration)
  max_cbs_iteration=as.integer(max_cbs_iteration)
  cv=as.integer(cv)
  cv_test_size=as.integer(floor(n/cv))
  
  result=get_h_candidate(X_training,cv,h_add,h_length)
  s=result$s
  X_training=X_training[s,]
  Y_training=Y_training[s,]
  h_candidate_matrix=result$h_candidate_matrix
  h_candidate_matrix[is.na(h_candidate_matrix)]=0
  
  result=.Fortran('CBS_density_simplex_LL_constraint',n=n,d_x_1=d,T=T,W_training=X_training,Y_training=Y_training,cv=cv,
                  cv_test_size=cv_test_size,add=add,Y_domain=Y_grid,I_or_F=as.integer(1),
                  max_smoothing_length=max(h_length),smoothing_length=h_length,smoothing_matrix=h_candidate_matrix,epsilon=epsilon,
                  max_cbs_iteration=max_cbs_iteration,max_sbf_iteration=max_iteration,
                  cbs_iteration=as.integer(0),sbf_iteration=array(as.integer(0),c(max_cbs_iteration,d,max(h_length),cv)),
                  optimal_smoothing=h_candidate_matrix[1,])
  
  return(list(optimal_h=result$optimal_smoothing,CBS_iteration=result$cbs_iteration,
              SBF_iteration=result$sbf_iteration))
}

# Function for fitting Bochner smooth backfitting
# X_test, X_training, Y_training are usual
# add: number of additional grid points for numerical integration
# Y_grid: common observed times for Y
# max_iteration: maximum iteration number for Bochner smooth bacfkfitting algorithm
# epsilon: convergence criterion for Bochner smooth backfitting algorithm
# h: ncol(X_training)-dimensional vector of bandwidths
SBF_density=function(X_test,X_training,Y_training,add=11,Y_grid=NULL,max_iteration=50,epsilon=10^-4,h)
{
  n=nrow(X_training)
  d=ncol(X_training)
  T=ncol(Y_training)
  if(is.matrix(X_test)) N=nrow(X_test)
  else N=as.integer(1)
  X_test=matrix(X_test,nrow=N,ncol=d)
  
  if(is.null(Y_grid)) Y_grid=as.double(seq(0,1,length=T))
  else Y_grid=as.double(Y_grid)
  add=as.integer(add)
  max_iteration=as.integer(max_iteration)
  h=as.double(h)
  
  g=integration.grid(X_test,add)
  ngrid=nrow(g)
  actual_iteration=as.integer(0)
  yhat=matrix(0,N,T)
  mhat=array(0,dim=c(N,d,T))
  
  result=.Fortran('SBF_density_simplex_LL_constraint',n0=N,n=n,d_x_1=d,T=T,W_test=X_test,W_training=X_training,Y_training=Y_training,
                  g_length=ngrid,g=g,Y_domain=Y_grid,I_or_F=as.integer(1),
                  smoothing_vector=h,epsilon=epsilon,max_sbf_iteration=max_iteration,sbf_iteration=actual_iteration,mhat=mhat,
                  yhat=yhat)
  
  return(list(Y_hat=result$yhat,m_hat=result$mhat,SBF_iteration=result$sbf_iteration))
  if(any(result$sbf_iteration==max_iteration)) print('Some smooth backfitting algorithm did not converge.')
}