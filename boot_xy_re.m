function [replication1,replication2] = boot_xy_re(r,i,k,t,total,yi_transpose,xi_transpose,tp1,constant,dof,bhat)
%Perform XY (pairs) bootstrap on clustered data

%r: number of boostrap replications
%i: number of cross section members
%k: number of independent variables in the model
%t: number of time periods
%total: i*t
%yi_transpose: i x t matrix of dependent variable data (row is cross
%section, column is time)
%xi_transpose: i x t matrix of independent variable data (row is cross
%section, column is time)
%tp1 = t+1
%constant: i*t vector of ones
%dof: n/(n-1) (used for crcme calculation
%bhat: estimate of slope coefficient from original data

%Returns 2xr matrix "replication1"
%replication1(1): ordered vector of bhat*s
%replication1(2): ordered vector of t statistics

%Returns 1xr matrix "replication2"
%replication2 : ordered vector of absolute value of t statistics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create storage for bootstrap replication
replication = zeros(2,r);

for bootrep = 1:r
    
    observations = randsample(i,i,true); %random sample of i cross section members with replacement
    data = [yi_transpose,xi_transpose]; %i x t matrix of observations
    sample = data(observations,:); %select all t observations with each cross section member
    yi_transpose_star = sample(:,1:t); %i x t matrix of dependent variable data
    xi_transpose_star = sample(:,tp1:end); %i x t matrix of independent variable data
    y_star = reshape(yi_transpose_star',total,1); %t*i x 1 vector of dependent variable data
    x_star = reshape(xi_transpose_star',total,1); %t*i x 1 vector of independent variable data (no constant)
    x1_star = [constant,x_star]; %Add constant to independent variable data
    xi_star = xi_transpose_star'; %t x i matrix of indepenent variable data without constant (used in crcme estimator)

    %Find estimates of Beta using OLS
    bhat_star = ((x1_star'*x1_star)^-1)*(x1_star'*y_star);
    residuals_boot = y_star - x1_star*bhat_star;
    residualsi_boot = reshape(residuals_boot,t,i); %t x i matrix of residuals used in crcme estimator
    
    %Calculate sum on t x k matrix of covariates for each cross section member
    %Used to calculate clustered robust standard errors
    sum = zeros(k,k); 
    for count = 1:i
        sum = sum + [ones(t,1),xi_star(:,count)]'*[ones(t,1),xi_star(:,count)];
    end
    x1iinv_star = sum^-1;
    
    cov_star = crcme(k,i,t,xi_star,residualsi_boot,dof,x1iinv_star);
    se_star = sqrt(cov_star);
    t_star = (bhat_star(k,1)-bhat(k,1))/se_star(k,k); %bootstrap t-statistic
    
    replication(1,bootrep) = bhat_star(k,1); %Store bhat*
    replication(2,bootrep) = t_star; %Store t-stat of bhat*

end

replication1 = sort(replication,2); %This is used for equal tailed percentile-t CIs

replication2 = abs(replication(2,:)); 
replication2 = sort(replication2,2); %This is used for symmetrical percentile-t CIs


end

