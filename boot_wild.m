function [replication1,replication2] = boot_wild(type_wild,wtest,f,r,n,residuals_h,x1bhat,x1inv,x1prime,x1,k,bias,ninv,bhat)
%Perform Wild Bootstrap
%type 1 = Rademacher
%type 2 = Mammen

%type_wild: 1 for Rademacher, 2 for Mammen
%wtest: .5 for Rademacher, (sqrt(5)+1)/(2*sqrt(5)); for Mammen
%f: [0,0] for Rademacher, [-(sqrt(5)-1)/2,(sqrt(5)+1)/2]; for Mammen
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

%Returns 2xr matrix "replication1"
%replication1(1): ordered vector of bhat*s
%replication1(2): ordered vector of t statistics

%Returns 1xr matrix "replication2"
%replication2 : ordered vector of absolute value of t statistics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create storage for bootstrap replication
replication = zeros(2,r);

%Create storage for perturbed residuals
residuals_star = zeros(n,1);


for bootrep = 1:r
    
    test = unifrnd(0,1,n,1);
    
    %Generate two-point distribution 
    if type_wild == 1
        for count = 1:n
            if (test(count,1) <= wtest)
                residuals_star(count,1) = -residuals_h(count,1);
            else
                residuals_star(count,1) = residuals_h(count,1);
            end
        end
        
    else
        for count = 1:n
        if (test(count,1) <= wtest)
            residuals_star(count,1) = f(1)*residuals_h(count,1);
        else 
            residuals_star(count,1) = f(2)*residuals_h(count,1);
        end
        end
        
    end
    
        
    
    y_star = x1bhat + residuals_star; %Calculate new dependent variable data y*  
    
    bhat_star = (x1inv)*(x1prime*y_star); %Calculate OLS estimate 
    residuals_boot = y_star - x1*bhat_star; %Calculate residuals
    
    %Compute robust standard error
    cov_star = hccme(k,n,bias,ninv,residuals_boot,x1inv,x1);
    se_star = sqrt(cov_star);
    t_star = (bhat_star(k,1)-bhat(k,1))/se_star(k,k); %bootstrap t-statistic
    
    replication(1,bootrep) = bhat_star(k,1); %Store bhat*
    replication(2,bootrep) = t_star; %Store t-stat of bhat*
   
end

replication1 = sort(replication,2); %This is used for equal tailed percentile-t CIs

replication2 = abs(replication(2,:)); 
replication2 = sort(replication2,2); %This is used for symmetrical percentile-t CIs


end