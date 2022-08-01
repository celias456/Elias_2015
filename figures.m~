
%% Load Data

clear all
clc;

load results_linear.csv 
data = results_linear;

linear_n = data(:,1);
linear_nominal = data(:,2);
linear_asymptotic = data(:,3);
linear_percentile_xy = data(:,4);
linear_percentile_t_xy = data(:,5);
linear_percentile_t2_xy = data(:,6);
linear_percentile_wild = data(:,7);
linear_percentile_t_wild = data(:,8);
linear_percentile_t2_wild = data(:,9);


load results_re.csv 
data = results_re;

re_i = data(:,1);
re_nominal = data(:,2);
re_asymptotic = data(:,3);
re_percentile_xy = data(:,4);
re_percentile_t_xy = data(:,5);
re_percentile_t2_xy = data(:,6);
re_percentile_wild = data(:,7);
re_percentile_t_wild = data(:,8);
re_percentile_t2_wild = data(:,9);


%% Generate Plots for Linear Model - XY Bootstrap

plot_linear_nominal = plot(linear_n,linear_nominal);
hold all
plot_linear_asymptotic = plot(linear_n,linear_asymptotic);
plot_linear_percentile = plot(linear_n,linear_percentile_xy);
plot_linear_percentile_t = plot(linear_n,linear_percentile_t_xy);
plot_linear_percentile_t2 = plot(linear_n,linear_percentile_t2_xy);
ylim([.45,.93])
xlim([10,300])
set(plot_linear_nominal,'Color','k','Linewidth',1.25);
set(plot_linear_asymptotic,'Color','k','LineStyle','--');
set(plot_linear_percentile,'Color',[.8 .8 .8],'LineStyle','-');
set(plot_linear_percentile_t,'Color','k','LineStyle','-');
set(plot_linear_percentile_t2,'Color','k','LineStyle',':');
xlabel('Sample Size');
ylabel('Empirical Coverage')
s1 = legend([plot_linear_asymptotic plot_linear_percentile plot_linear_percentile_t plot_linear_percentile_t2],'Asymptotic','Percentile','Equal-Tailed Percentile-t','Symmetric Percentile-t','location','SouthEast');
hold off


%% Generate Plots for Linear Model - Wild Bootstrap

plot_linear_nominal = plot(linear_n,linear_nominal);
hold all
plot_linear_asymptotic = plot(linear_n,linear_asymptotic);
plot_linear_percentile = plot(linear_n,linear_percentile_wild);
plot_linear_percentile_t = plot(linear_n,linear_percentile_t_wild);
plot_linear_percentile_t2 = plot(linear_n,linear_percentile_t2_wild);
ylim([.45,.93])
xlim([10,300])
set(plot_linear_nominal,'Color','k','Linewidth',1.25);
set(plot_linear_asymptotic,'Color','k','LineStyle','--');
set(plot_linear_percentile,'Color',[.8 .8 .8],'LineStyle','-');
set(plot_linear_percentile_t,'Color','k','LineStyle','-');
set(plot_linear_percentile_t2,'Color','k','LineStyle',':');
xlabel('Sample Size');
ylabel('Empirical Coverage')
s2 = legend([plot_linear_asymptotic plot_linear_percentile plot_linear_percentile_t plot_linear_percentile_t2],'Asymptotic','Percentile','Equal-Tailed Percentile-t','Symmetric Percentile-t','location','SouthEast');
hold off



%% Generate Plots for RE Model - XY Bootstrap

plot_re_nominal = plot(re_i,re_nominal);
hold all
plot_re_asymptotic = plot(re_i,re_asymptotic);
plot_re_percentile = plot(re_i,re_percentile_xy);
plot_re_percentile_t = plot(re_i,re_percentile_t_xy);
plot_re_percentile_t2 = plot(re_i,re_percentile_t2_xy);
ylim([.65,.93])
xlim([5,60])
set(plot_re_nominal,'Color','k','Linewidth',1.25);
set(plot_re_asymptotic,'Color','k','LineStyle','--');
set(plot_re_percentile,'Color',[.8 .8 .8],'LineStyle','-');
set(plot_re_percentile_t,'Color','k','LineStyle','-');
set(plot_re_percentile_t2,'Color','k','LineStyle',':');
xlabel('Cross Section Size');
ylabel('Empirical Coverage')
r1 = legend([plot_re_asymptotic plot_re_percentile plot_re_percentile_t plot_re_percentile_t2],'Asymptotic','Percentile','Equal-Tailed Percentile-t','Symmetric Percentile-t','location','SouthEast');
hold off


%% Generate Plots for RE Model - Wild Bootstrap



plot_re_nominal = plot(re_i,re_nominal);
hold all
plot_re_asymptotic = plot(re_i,re_asymptotic);
plot_re_percentile = plot(re_i,re_percentile_wild);
plot_re_percentile_t = plot(re_i,re_percentile_t_wild);
plot_re_percentile_t2 = plot(re_i,re_percentile_t2_wild);
ylim([.65,.93])
xlim([5,60])
set(plot_re_nominal,'Color','k','Linewidth',1.25);
set(plot_re_asymptotic,'Color','k','LineStyle','--');
set(plot_re_percentile,'Color',[.8 .8 .8],'LineStyle','-');
set(plot_re_percentile_t,'Color','k','LineStyle','-');
set(plot_re_percentile_t2,'Color','k','LineStyle',':');
xlabel('Cross Section Size');
ylabel('Empirical Coverage')
r2 = legend([plot_re_asymptotic plot_re_percentile plot_re_percentile_t plot_re_percentile_t2],'Asymptotic','Percentile','Equal-Tailed Percentile-t','Symmetric Percentile-t','location','SouthEast');
hold off

