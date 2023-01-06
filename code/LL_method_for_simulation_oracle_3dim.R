# LL method

# Load necessary package
library(pracma)

# Load necessary source file
source('C:/Users/u0124228/Desktop/Library/Research/Smooth Backfitting for Local Linear Hilbertian Additive Regression/R_functions_for_LL_method_for_density_responses_constraint_weight.R')
load(file='C:/Users/u0124228/Desktop/Library/Research/Smooth Backfitting for Local Linear Hilbertian Additive Regression/Oracle_data_3dim.Rdata')

# Load necessary dll files
dyn.load('C:/Users/u0124228/Desktop/Final Fortran Codes/New Codes for Local Linear/CBS_density_simplex_LL_constraint.dll')
dyn.load('C:/Users/u0124228/Desktop/Final Fortran Codes/New Codes for Local Linear/SBF_density_simplex_LL_constraint_weight.dll')

# Function for calculating the norm of density space
density_norm=function(d1,d2,time_vector)
{
  T=length(time_vector)
  first_integral=c()
  for(t in 1:T)
  {
    first_integral[t]=trapz(time_vector,(log(d1/d1[t])-log(d2/d2[t]))^2)
  }
  trapz(time_vector,first_integral)/2
}

# Defining component maps
new_f1=function(x,t)
{
  exp(-(x-1/2)*t)
}

new_f2=function(x,t)
{
  exp(-2*(x^2-1/3)*t^2)
}

new_f3=function(x,t)
{
  exp(-3*(x^3-1/4)*t^3)
}

# Function for results
# X: array for predictor of dimension (number of monte-carlo simulation(M) by training sample size(n) by dimension of predictor(d))
# Y: array for response of dimension (number of monte-carlo(M) by training sample size(n) by number of evaluation time points(T))
# time_vector: vector of evaluation time points
# h_add and h_length: d-dimensional vectors 
# seq(min_h[j],min_h[j]+h_add[j],length=h_length[j]) will be the candidate set of j-th bandwidth
# for some small bandwidth min_h[j]
LL=function(X,Y,time_vector,h_length,h_add)
{
  M=dim(X)[1]
  n=dim(X)[2]
  d=dim(X)[3]
  T=length(time_vector)
  
  # Integration grid of x for IMSE
  x_add=51
  grid=seq(0,1,length=x_add)
  grid_matrix=matrix(rep(grid,d),ncol=d)
  
  # m_true
  m_true=array(,dim=c(d,x_add,T))
  for(t in 1:T)
  {
    m_true[1,,t]=new_f1(grid,time_vector[t])
    m_true[2,,t]=new_f2(grid,time_vector[t])
    m_true[3,,t]=new_f3(grid,time_vector[t])
  }
  for(l in 1:x_add)
  {
    m_true[1,l,]=m_true[1,l,]/trapz(time_vector,m_true[1,l,])
    m_true[2,l,]=m_true[2,l,]/trapz(time_vector,m_true[2,l,])
    m_true[3,l,]=m_true[3,l,]/trapz(time_vector,m_true[3,l,])
  }
  
  # Compute m_hat
  new_weight=array(,dim=c(n,x_add,d))
  before_sum=array(,dim=c(n,x_add,d,T))
  m_hat=array(,dim=c(M,d,x_add,T))
  for(j in 1:M)
  {
    print(j)
    optimal_h=CBS_density(X[j,,],Y[j,,],add=2,Y_grid=time_vector,h_length=h_length,h_add=h_add)$optimal_h
    final_weight=SBF_density(grid_matrix,X[j,,],Y[j,,],add=0,h=optimal_h,Y_grid=time_vector)$final_weight
    for(i in 1:n)
    {
      for(k in 1:d)
      {
        new_weight[i,,k]=final_weight[i,,k]-trapz(grid,final_weight[i,,k])
        for(l in 1:x_add)
        {
          before_sum[i,l,k,]=Y[j,i,]^{new_weight[i,l,k]/n}
        }
      }
    }
    for(k in 1:d)
    {
      for(l in 1:x_add)
      {
        prods=colProds(before_sum[,l,k,])
        m_hat[j,k,l,]=prods/trapz(time_vector,prods)
      }
    }
  }
  
  # Compute average of m_hat
  avg_m_hat=array(,dim=c(d,x_add,T))
  for(k in 1:d)
  {
    for(l in 1:x_add)
    {
      prods=colProds(m_hat[,k,l,]^{1/M})
      avg_m_hat[k,l,]=prods/trapz(time_vector,prods)
    }
  }
  
  bias_vector=matrix(0,d,x_add)
  variance_vector=matrix(0,d,x_add)
  bias=c()
  variance=c()
  norm_vector=c()
  norm_temporary=c()
  
  # ISB
  for(k in 1:d)
  {
    for(l in 1:x_add)
    {
      bias_vector[k,l]=density_norm(m_true[k,l,],avg_m_hat[k,l,],time_vector)
    }
    bias[k]=trapz(grid,bias_vector[k,])
  }
  
  # IV
  for(k in 1:d)
  {
    for(l in 1:x_add)
    {
      for(j in 1:M)
      {
        norm_temporary[j]=density_norm(m_hat[j,k,l,],avg_m_hat[k,l,],time_vector)
      }
      variance_vector[l]=mean(norm_temporary)
    }
    variance[k]=trapz(grid,variance_vector)
  }
  
  return(list(bias=bias,variance=variance,bias_vector=bias_vector,variance_vector=variance_vector,m_hat=m_hat,new_weight=new_weight,h=optimal_h))
}

d=dim(data$X)[3]
result=LL(data$X,data$Y,seq(-0.5,0.5,length=101),h_length=rep(21,d),h_add=rep(0.2,d))