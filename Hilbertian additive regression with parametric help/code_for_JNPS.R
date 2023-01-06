library(ks)
library(pracma)

# Load necessary source file, the collection of R functions that are used to compute the local linear SBF estimators for density-responses
source('D://R_functions_for_LL_method_for_density_responses_constraint.R')

# Fortran dll files performing the SBF algorithm
dyn.load('D://SBF_density_simplex_LL_constraint_derivative.dll') 
dyn.load('D://SBF_density_simplex_LL_constraint.dll') 

# Fortran dll file performing the CBS bandwidth selection for the local linear SBF method
dyn.load('D://CBS_density_simplex_LL_constraint.dll')


# Function for computing the local linear SBF estimators of mean regression map and its component maps
SBF_density_derivative=function(X_test,X_training,Y_training,add=11,Y_grid=NULL,max_iteration=50,epsilon=10^-4,h)
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
  mdhat=array(0,dim=c(N,d,T))
  
  result=.Fortran('SBF_density_simplex_LL_constraint_derivative',n0=N,n=n,d_x_1=d,T=T,W_test=X_test,W_training=X_training,Y_training=Y_training,
                  g_length=ngrid,g=g,Y_domain=Y_grid,I_or_F=as.integer(1),
                  smoothing_vector=h,epsilon=epsilon,max_sbf_iteration=max_iteration,sbf_iteration=actual_iteration,mhat=mhat,
                  yhat=yhat,mdhat=mdhat)
  
  return(list(Y_hat=result$yhat,m_hat=result$mhat,m_d_hat=result$mdhat,SBF_iteration=result$sbf_iteration))
  if(any(result$sbf_iteration==max_iteration)) print('Some smooth backfitting algorithm did not converge.')
}

# Function for calculating the distance between two densities
density_distance=function(d1,d2,time_vector)
{
  T=length(time_vector)
  first_integral=c()
  for(t in 1:T)
  {
    first_integral[t]=trapz(time_vector,(log(d1/d1[t])-log(d2/d2[t]))^2)
  }
  trapz(time_vector,first_integral)/2
}

# Function for calculating the inner product of two densities
density_inprod=function(d1,d2,time_vector)
{
  T=length(time_vector)
  temp1=c()
  for(t in 1:T)
  {
    temp1[t]=trapz(time_vector,(log(d1/d1[t])*log(d2/d2[t])))
  }
  trapz(time_vector,temp1)/2
}

# Data pre-processing: generating smooth densities from raw data
library(accelmissing)
data(acceldata2)
value=as.matrix(acceldata2$PA)

length(which(value>5000))/length(which(value>0))

adult=which(acceldata2$demo[,2]>=18)

n=184
tresh=5000
eval_number=51
eval_points=seq(0,log(tresh),length=eval_number)
N=rep(0,n)
Y=matrix(0,nrow=n,ncol=eval_number)
for(i in adult)
{
  x=as.vector(value[((i-1)*7+1):(i*7),])
  x=x[which(x<=tresh)]
  x=log(x[which(x>0)])
  N[i]=length(x)
  bw=hscv(x)
  Y_temporary=kde(x,h=bw,eval.points=eval_points)$estimate
  if(length(which(Y_temporary<=0))>0) Y_temporary[which(Y_temporary<=0)]=min(Y_temporary[which(Y_temporary>0)])
  Y[i,]=Y_temporary/trapz(eval_points,Y_temporary)
}

Y=Y[which(rowSums(Y)>0),]
n=nrow(Y)

X1=log(acceldata2$demo[adult,2])
X2=log(acceldata2$demo[adult,4])
X1=(X1-min(X1))/(max(X1)-min(X1))
X2=(X2-min(X2))/(max(X2)-min(X2))
X=cbind(X1,X2)
d=ncol(X)

############ Computation of Hilbertian Local Linear SBF ####################

optimal_h_LL=matrix(0,nrow=n,ncol=d)
Y_hat_HLL=matrix(0,nrow=n,ncol=eval_number)
Y_mhat_HLL=array(0,c(d,n,eval_number))
error_LL=c()

# Get ASPE for Hilbertian LL-SBF
for(k in 1:n)
{
  print(k)
  X_test=X[k,]
  Y_test=Y[k,]
  X_training=X[-k,]
  Y_training=Y[-k,]
  optimal_h_LL[k,]=CBS_density(X_training,Y_training,h_add=c(0.4,0.4),h_length=c(41,41),cv=n-1,Y_grid=eval_points)$optimal_h
  Y_hat=SBF_density(X_test,X_training,Y_training,h=optimal_h_LL[k,],Y_grid=eval_points)$Y_hat
  Y_hat_HLL[k,]<-c(Y_hat)
  
  error_LL[k]=density_distance(Y_test,Y_hat,eval_points)
  print(mean(error_LL))
}

mean(error_LL)/log(tresh) # ASPE for Hilbertian LL-SBF


############## For Hilbertian LL SBF with Para-Help ####################
## compute Hilbertian hat m_j(Xj)(cdot) estimates
# choosing gj=hat m_{Y,j} : the estimated components in the additive decomposition of local linear SBF est.

Y_mhat_HLL=array(0,c(d,n,eval_number)) # leave-one-out version hat mg_j
for(k in 1:n)
{
  print(k)
  X_test=X[k,]
  Y_test=Y[k,]
  X_training=X[-k,]
  Y_training=Y[-k,]
 
  X1_test=cbind(c(rep(X[k,1],n-1)),X[-k,2])
  result_X1=SBF_density_derivative(matrix(X1_test,(n-1),d),X_training,Y_training,Y_grid=eval_points,h=optimal_h_LL[k,])
  g1_x1<-apply(result_X1$m_hat[,1,],2,mean)
  X2_test=cbind(X[-k,1],c(rep(X[k,2],n-1)))
  result_X2=SBF_density_derivative(matrix(X2_test,(n-1),d),X_training,Y_training,Y_grid=eval_points,h=optimal_h_LL[k,])
  g2_x2<-apply(result_X2$m_hat[,2,],2,mean)
  
  Y_mhat_HLL[1,k,]<- g1_x1/trapz(eval_points,g1_x1)
  Y_mhat_HLL[2,k,]<- g2_x2/trapz(eval_points, g2_x2)
} # end of getting g_j 


## gj(Xij) = hat m_j(Xij)
ft_g<-Y_mhat_HLL # new responses for G1(X1) and G2(X2) leave-one-out version

## compute Hilbertian local linear SBF estimate hat m_{gj(Xj)}
## m_hat_G : hat m_{gj(Xj)} Hilbertian local linear SBF with {(Xi,gj(Xij))} 
## compute error(g(X)) and error(Y)
optimal_g1h_LL=matrix(0,nrow=n,ncol=d)
optimal_g2h_LL=matrix(0,nrow=n,ncol=d)
m_hat_G=array(0,c(d,n,eval_number))
error_gX<-array(0,c(d,n,eval_number))
error_Y<-matrix(0,nrow=n,ncol=eval_number)
for(k in 1:n)
{
  print(k)
  X_test=X[k,]
  G_test=ft_g[,k,]
  X_training=X[-k,]
  G_training=ft_g[,-k,]
  optimal_g1h_LL[k,]=CBS_density(X_training,G_training[1,,],h_add=c(0.4,0.4),h_length=c(41,41),cv=n-1,Y_grid=eval_points)$optimal_h
  m_hat_g1=SBF_density(X_test,X_training,G_training[1,,],h=optimal_g1h_LL[k,],Y_grid=eval_points)$Y_hat
  m_hat_G[1,k,]<-c(m_hat_g1)#/trapz(eval_points,m_hat_g1)
  optimal_g2h_LL[k,]=CBS_density(X_training,G_training[2,,],h_add=c(0.4,0.4),h_length=c(41,41),cv=n-1,Y_grid=eval_points)$optimal_h
  m_hat_g2=SBF_density(X_test,X_training,G_training[2,,],h=optimal_g2h_LL[k,],Y_grid=eval_points)$Y_hat
  m_hat_G[2,k,]<-c(m_hat_g2)#/trapz(eval_points,m_hat_g2)

  # error 
  error_gX[1,k,]<-(ft_g[1,k,]/m_hat_g1)/trapz(eval_points,(ft_g[1,k,]/m_hat_g1))
  error_gX[2,k,]<-(ft_g[2,k,]/m_hat_g2)/trapz(eval_points,(ft_g[2,k,]/m_hat_g2))

  error_Y[k,]<-(Y[k,]/Y_hat_HLL[k,])/trapz(eval_points,(Y[k,]/Y_hat_HLL[k,]))
  }
 

# compute theta according to eq(1.1);-eq(2.7) revision on page 8

t11=sum(sapply(1:n, function(k) density_inprod(error_gX[1,k,],error_gX[1,k,],eval_points)))
t12=sum(sapply(1:n, function(k) density_inprod(error_gX[1,k,],error_gX[2,k,],eval_points)))
t22=sum(sapply(1:n, function(k) density_inprod(error_gX[2,k,],error_gX[2,k,],eval_points)))

th_1=sum(sapply(1:n, function(k) density_inprod(error_gX[1,k,],error_Y[k,],eval_points)))
th_2=sum(sapply(1:n, function(k) density_inprod(error_gX[2,k,],error_Y[k,],eval_points)))

theta<-solve(matrix(c(t11,t12,t12,t22),2,2))%*%c(th_1,th_2)



# Get ASPE for Hilbertian LL-SBF with para-help
Y_hat_HPH<-matrix(0,nrow=n,ncol=eval_number)
error_PH<-c()

for(k in 1:n){
  temp1<-error_gX[1,k,]^theta[1]/trapz(eval_points,(error_gX[1,k,])^theta[1])
  temp2<-error_gX[2,k,]^theta[2]/trapz(eval_points,(error_gX[2,k,])^theta[2])
  tempk<-Y_hat_HLL[k,]*temp1*temp2
  Y_hat_HPH[k,]<-tempk/trapz(eval_points,tempk)
  error_PH[k]<-density_distance(Y[k,],Y_hat_HPH[k,],eval_points)
  print(mean(error_PH))
}

mean(error_PH)/log(tresh)

mean(error_LL)/log(tresh)

