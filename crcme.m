function [cov] = crcme(k,i,t,xi,residualsi,dof,x1iinv)
%Clustered Robust Covariance Matrix Estimator

%k: number of independent variables
%i: number of cross section members
%xi: t x i matrix of covariates without constant (row is time, column is cross section member) 
%residualsi: t x i matrix of residuals (row is time, column is cross section
%member) 
%dof: degrees of freedom correction
%x1iinv: outside matrix sum used in crcme (See Greene(2007))

sum = zeros(k,k);
for count = 1:i
    sum = sum + [ones(t,1),xi(:,count)]'*residualsi(:,count)*residualsi(:,count)'*[ones(t,1),xi(:,count)];
end
cov = dof.*x1iinv*sum*x1iinv;



end

