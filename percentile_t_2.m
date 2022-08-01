function [ci] = percentile_t_2(replication,quantile,bhat_se,k,bhat,Beta)
%Percentile-t Confidence Interval - Symmetric

%replication: ordered vector of absolute value of bootstrapped t statistics
%quantile: appropriate quantile for CI construction
%bhat_se: standard error of bhat
%k: number of independent variables
%bhat: estimate of slope coefficient from original data
%Beta: True value of the slope coefficient

%Returns 3x1 vector 'ci'
%ci(1): lower bound of confidence interval
%ci(2): upper bound of confidence interval
%ci(3): equals 1 if true slope is in inverval, 0 otherwise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ci = zeros(3,1);
cv = replication(quantile);
se = bhat_se(k,k);

%Form percentile t confidence interval
ci(1) = bhat(k,1) - cv*se; %lower bound
ci(2) = bhat(k,1) + cv*se; %upper bound

%Determine if percentile t confidence interval contains true Beta 1
if ((Beta(k,1) >= ci(1)) && (Beta(k,1) <= ci(2)))
    ci(3) = 1;


end

