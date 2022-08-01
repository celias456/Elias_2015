function [replication] = boot_residual(r,n,residuals_h,x1bhat,x1inv,x1prime,x1,k,bias,ninv,bhat)
%Perform Nonparametric Residual Bootstrap

%r: number of bootstrap replications
%n: sample size
%residuals_h: residuals from original sample removed of small sample bias
%and heteroscedasticity
%x1bhat: x1*bhat
%x1inv: (x1'*x1)^-1
%x1prime: x1'
%x1: independent variable data with constant
%k: number of independent variables
%bias: n/(n-k)
%ninv: 1/n
%bhat: ols estimate of Beta from original data

%Returns 2xr matrix
%replication(1): ordered vector of bhat*s
%replication(2): ordered vector of t statistics


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create storage for bootstrap replication
replication = zeros(2,r);


for bootrep = 1:r
    
    sample = randsample(n,n,true); %random sample of n observations with replacement
    residuals_star = residuals_h(sample); %Take random sample with replacement from LS residuals corrected for small sample bias and heteroscedasticity
    y_star = x1bhat + residuals_star; %Calculate new dependent variable data Ystar
    bhat_star = (x1inv)*(x1prime*y_star); %Calculate OLS estimate with Ystar data
    residuals_boot = y_star - x1*bhat_star; %Calculate residuals
    
    %Compute robust standard error
    cov_star = hccme(k,n,bias,ninv,residuals_boot,x1inv,x1);
    se_star = sqrt(cov_star);
    t_star = (bhat_star(k,1)-bhat(k,1))/se_star(k,k); %bootstrap t-statistic
    
    replication(1,bootrep) = bhat_star(k,1); %Store bhat*
    replication(2,bootrep) = t_star; %Store t-stat of bhat*
    
    
end

replication = sort(replication,2);

end

