function [ci] = percentile_t(bhat,k,replication,low,up,bhat_se,Beta)
%Percentile-t Confidence Interval - Equal Tailed

%bhat: estimate of slope coefficient from original data
%k: number of independent variables
%replication: ordered vector of bootstrapped t statistics
    %replication(1,:): ordered vector of bootstrapped beta-hats
    %replication(2,:): ordered vector of bootstrapped t statistics
%low: Quantile used for upper bound of confidence interval
%up: Quantile used for lower bound of confidence interval
%bhat_se: standard error of bhat
%Beta: True value of slope coefficient

%Returns 3x1 vector 'percentile_t'
%percentile_t(1): lower bound of confidence interval
%percentile_t(2): upper bound of confidence interval
%percentile_t(3): equals 1 if true slope is in inverval, 0 otherwise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ci = zeros(3,1);

%Form percentile t confidence interval
ci(1) = bhat(k,1) - (replication(2,up)*(bhat_se(k,k))); %lower bound
ci(2) = bhat(k,1) - (replication(2,low)*(bhat_se(k,k))); %upper bound

%Determine if percentile t confidence interval contains true Beta 1
if ((Beta(k,1) >= ci(1)) && (Beta(k,1) <= ci(2)))
    ci(3) = 1;

end

