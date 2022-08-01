%Random Effects Model

clear all
clc;

%% Set Design Parameters and Preliminaries
i = 26; %Set value for cross section 
t = 5; %Set value for number of time periods
total = i*t; %Number of total observations to generate
criticalvalue = 1.65685; %Critical value for confidence interval
type_error_dist = 1; %1 for standard normal error, 2 for student's t error
type_error_variance = 2; %1 for homoscedasticity, %2 for heteroscedasticity
type_wild = 1; %1 = Rademacher, 2 = Mammen
seedstate = 7; %Set seed state

r = 999; %Number of bootstrap repititions
m = 1000; %Number of monte carlo repititions
dof = i/(i-1);
tp1 = t+1;

stream = RandStream('mt19937ar','Seed',seedstate);
RandStream.setDefaultStream(stream);

Beta = [1;1];  %True value of beta
k = length(Beta); %Determine the number of independent variables in the model

%Confidence level 
alpha = .1;
confidence = 1-alpha;
low = (alpha/2)*(r+1);
up = (1-(alpha/2))*(r+1);
quantile = (1-alpha)*(r+1);

%Wild bootstrap characteristics
if type_wild == 1
    wtest = .5;
    f = [0,0];
else
    wtest = (sqrt(5)+1)/(2*sqrt(5));
    f = [-(sqrt(5)-1)/2,(sqrt(5)+1)/2];
end

%% Independent variable data
load xdata.csv 
constant = ones(total,1);
x = xdata(1:total,2);
xi = reshape(x,t,i); %t x i matrix; each column is cross section member, row is time
xi_transpose = xi'; %i x t matrix; row is cross section, column is time
x1 = xdata(1:total,:); % x data with constant
x1inv = (x1'*x1)^-1;
x1prime = x1';
x1Beta = x1*Beta;
person = sort(repmat((1:1:i)',t,1)); %cross section index
time = repmat((1:1:t)',i,1); %time index
data = [person,time,x1];

% Perform calculations on X data to remove heteroscedasticity and small sample 
% bias
h = x1*(x1inv)*x1';
hi = zeros(total,1);
for count_h = 1:total
    hi(count_h,1) = sqrt(1-h(count_h,count_h));
end

%Calculate sum of t x k matrix of covariates for each cross section member
%Used to calculate clustered robust standard errors
sumxi = zeros(k,k); 
for count_xi = 1:i
    sumxi = sumxi + [ones(t,1),xi(:,count_xi)]'*[ones(t,1),xi(:,count_xi)];
end
x1iinv = sumxi^-1;

%% Error term
if type_error_dist == 1 %Gaussian error term
    errormean = zeros(1,t); %Mean of the error term
else %Student's t error term
    degreesoffreedom = 3;
end

muvariance = 1; %Variance of the common error term
epsilonvariance = 9; %Variance of the total error term
capsigma = epsilonvariance.*eye(t)+muvariance.*(ones(t,1)*ones(1,t)); %Form cov matrix for cross section members

if type_error_variance == 1 %Homoscedasticity
    sigma = ones(total,1);
else %Heteroscedasticity
    sigma = x;
end

%% Initialize variables for simulation

%t x i matrix for error term storage
ei = zeros(t,i); 

%Storage for confidence interval bounds (%first row is lower bound, second
%row is upper bound, third row records if Beta is in ci (1=yes,0=no))
ci_asy = zeros(3,m); 

ci_xy_percentile = zeros(3,m);
ci_xy_percentile_t = zeros(3,m);
ci_xy_percentile_t_2 = zeros(3,m);

ci_wild_percentile = zeros(3,m);
ci_wild_percentile_t = zeros(3,m);
ci_wild_percentile_t_2 = zeros(3,m);

%% Start Monte Carlo Repetitions

for montecarlocount = 1:m;

%% Generate error term and dependent variable data

if type_error_dist == 1 %Gaussian error
    for count_error = 1:i
        ei(:,count_error) = mvnrnd(errormean,capsigma)'; %t x i matrix of error terms
    end
else %Student's t error
    for count_error = 1:i
        ei(:,count_error) = mvtrnd(capsigma,degreesoffreedom)'; %t x i matrix of error terms
    end
end
    
e = reshape(ei,total,1); %i*t x 1 matrix of error terms
    
y = x1Beta + sigma.*e; %Calculate dependent variable Y  
yi = reshape(y,t,i); %t x i matrix; each column is cross section member, row is time
yi_transpose = yi'; %i x t matrix; row is cross section, column is time
    
%% Perform ols regression
bhat = x1inv*(x1prime*y); %Calculate LS estimate
x1bhat = x1*bhat;
residuals = y - x1bhat; %Calculate residuals
residuals_h = residuals./hi; %Remove small sample bias and heteroscedasticity from residuals
residualsi = reshape(residuals,t,i); %t x i matrix; each column is cross section member, row is time
residualsi_h = reshape(residuals_h,t,i); %t x i matrix; each column is cross section member, row is time
residualsi_h_transpose = residualsi_h'; %i x t matrix; row is cross section, column is time

%Calculate robust standard errors
bhat_cov = crcme(k,i,t,xi,residualsi,dof,x1iinv);
bhat_se = sqrt(bhat_cov);

%% Generate asymptotic confidence interval of Beta slope coefficient

ci_asy(1,montecarlocount) = bhat(k,1) - criticalvalue*bhat_se(k,k); %Beta 1 lower bound
ci_asy(2,montecarlocount) = bhat(k,1) + criticalvalue*bhat_se(k,k); %Beta 1 upper bound

%Determine if Asymptotic confidence interval contains true Beta 1
if ((Beta(k,1) >= ci_asy(1,montecarlocount)) && (Beta(k,1) <= ci_asy(2,montecarlocount)))
    ci_asy(3,montecarlocount) = 1;
end

%% XY Bootstrap

[replication_xy,replication2_xy] = boot_xy_re(r,i,k,t,total,yi_transpose,xi_transpose,tp1,constant,dof,bhat);

ci_xy_percentile(:,montecarlocount) = percentile(replication_xy,low,up,Beta,k);
ci_xy_percentile_t(:,montecarlocount) = percentile_t(bhat,k,replication_xy,low,up,bhat_se,Beta);
ci_xy_percentile_t_2(:,montecarlocount) = percentile_t_2(replication2_xy,quantile,bhat_se,k,bhat,Beta);


%% Wild Bootstrap

[replication_wild,replication2_wild] = boot_wild_re(type_wild,wtest,f,r,i,k,t,total,residualsi_h,x1bhat,x1inv,x1prime,x1,xi,dof,x1iinv,bhat);

ci_wild_percentile(:,montecarlocount) = percentile(replication_wild,low,up,Beta,k);
ci_wild_percentile_t(:,montecarlocount) = percentile_t(bhat,k,replication_wild,low,up,bhat_se,Beta);
ci_wild_percentile_t_2(:,montecarlocount) = percentile_t_2(replication2_wild,quantile,bhat_se,k,bhat,Beta);
    
end

%% Find Empirical Coverage Rates

count_ci_asy = sum(ci_asy(3,:));

count_ci_xy_percentile = sum(ci_xy_percentile(3,:));
count_ci_xy_percentile_t = sum(ci_xy_percentile_t(3,:));
count_ci_xy_percentile_t_2 = sum(ci_xy_percentile_t_2(3,:));

count_ci_wild_percentile = sum(ci_wild_percentile(3,:));
count_ci_wild_percentile_t = sum(ci_wild_percentile_t(3,:));
count_ci_wild_percentile_t_2 = sum(ci_wild_percentile_t_2(3,:));

ec_asy = count_ci_asy/m;

ec_xy_percentile = count_ci_xy_percentile/m;
ec_xy_percentile_t = count_ci_xy_percentile_t/m;
ec_xy_percentile_t_2 = count_ci_xy_percentile_t_2/m;

ec_wild_percentile = count_ci_wild_percentile/m;
ec_wild_percentile_t = count_ci_wild_percentile_t/m;
ec_wild_percentile_t_2 = count_ci_wild_percentile_t_2/m;


empiricalcoverage = struct('i',i,'sample_size',total,'critical_value',criticalvalue,'asymptotic',ec_asy,'xy_percentile',ec_xy_percentile,'xy_percentile_t',ec_xy_percentile_t,'xy_percentile_t_2',ec_xy_percentile_t_2,'wild_percentile',ec_wild_percentile,'wild_percentile_t',ec_wild_percentile_t,'wild_percentile_t_2',ec_wild_percentile_t_2);
