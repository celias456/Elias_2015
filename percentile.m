function [ci] = percentile(replication,low,up,Beta,k)
%Percentile Confidence Interval

%replication: 2 x r matrix
    %replication(1,:): ordered vector of bootstrapped beta-hats
    %replication(2,:): ordered vector of bootstrapped t statistics
%low: Quantile used for lower bound of confidence interval
%up: Quantile used for upper bound of confidence interval
%Beta: True value of slope coefficient
%k: Number of independent variables

%Returns 3x1 vector 'ci'
%ci(1): lower bound of confidence interval
%ci(2): upper bound of confidence interval
%ci(3): equals 1 if true slope coefficient is in interval, 0 otherwise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ci = zeros(3,1);

%Form percentile confidence interval
ci(1) = replication(1,low); %lower bound
ci(2) = replication(1,up); %upper bound

%Determine if percentile confidence interval contains true Beta 1
if ((Beta(k,1) >= ci(1)) && (Beta(k,1) <= ci(2)))
    ci(3) = 1;
end


end