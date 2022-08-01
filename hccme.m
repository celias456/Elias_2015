function [cov] = hccme(k,n,bias,ninv,residuals,x1inv,x1)
%Heteroscedastic consistent covariance matrix estimator
%Calculates a HCCME for given data

%k: number of independent variables
%n: sample size
%bias: bias correction (n/(n-k))
%ninv: 1/n
%residuals: ols residuals
%x1inv: (x1'x1)^-1
%x1: independent variable data with constant

%Returns k x k hccme covariance matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sum = zeros(k,k);

for count = 1:n
    sum = sum + (residuals(count)^2)*((x1(count,:)')*x1(count,:));
end
sum = ninv*sum;

cov = (bias)*(n)*(x1inv)*sum*(x1inv);

end

