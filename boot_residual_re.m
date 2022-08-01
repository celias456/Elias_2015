function [replication] = boot_residual_re(r,i,k,t,total,residualsi_h_transpose,x1bhat,x1inv,x1prime,x1,xi,dof,x1iinv,bhat)
%Performs Nonparametric Residual Bootstrap on Clustered Data

%r: number of boostrap replications
%i: number of cross section members
%k: number of independent variables in the model
%t: number of time periods
%total: i*t
%residuals_h_transpose: i x t matrix of residuals corrected for
%heteroscedasticity and small sample bias (row is cross section, column is
%time)
%x1bhat: x1*bhat
%x1inv: (x1'*x1)^-1
%x1prime: x1'
%x1: independent variable data with constant
%xi: t x i matrix of independent variable data (row is time, column is
%cross section member)
%dof: n/(n-1) (used for crcme calculation
%x1iinv: outer matrix used in crcme calculation
%bhat: estimate of slope coefficient from original data


%Returns 2xr matrix
%replication(1): ordered vector of bhat*s
%replication(2): ordered vector of t statistics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create storage for bootstrap replication
replication = zeros(2,r);

for bootrep = 1:r
    sample = randsample(i,i,true); %Random sample of i cross section members with replacement
    residualsi_transpose_star = residualsi_h_transpose(sample,:); %Take sample of residuals with replacement
    residuals_star = reshape(residualsi_transpose_star',total,1); %t*i vector of sampled residuals
    y_star = x1bhat+residuals_star; %dependent variable data

    %Find estimates of Beta using OLS
    bhat_star = (x1inv)*(x1prime*y_star);
    residuals_boot = y_star - x1*bhat_star; %t*i matrix of residuals
    residualsi_boot = reshape(residuals_boot,t,i); %t x i matrix of residuals used in crcme estimator

    cov_star = crcme(k,i,t,xi,residualsi_boot,dof,x1iinv);
    se_star = sqrt(cov_star);
    t_star = (bhat_star(k,1)-bhat(k,1))/se_star(k,k); %bootstrap t-statistic
    
    replication(1,bootrep) = bhat_star(k,1); %Store bhat*
    replication(2,bootrep) = t_star; %Store t-stat of bhat*  
end

replication = sort(replication,2);

end

