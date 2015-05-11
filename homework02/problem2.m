%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
% problem2.m                                                                   %
%                                                                              %
% These codes illustrate an implementation of fitting a curve to data using    %
% linear regression by updating the canonical information for linear regression%
% presented in Lectures 04a - 04b of Math 491 - Spring 2015.                   %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
F = @(x)( [ones(size(x)), x] ); % Simple linear regression functions f1(x) = 1, f2(x) = x
a = [2 3]'; % Parameters to fit
sigma2 = .2; % Epsilon variance

% Something more complicated is possible
%F = @(x)( [ones(size(x)), x, x.^2, sin(x)] ); 
%a = [3 1 4 1]'; 

n1 = 10; % How many points to fit initially
n2 = 20; % How many additional points at to fit at second stage
n3 = 50; % How many additional points at to fit at third stage

x1 = rand(n1,1);% Generate n1 points in [0,1]
y1 = F(x1)*a + sqrt(sigma2)*randn(n1,1);% Generate noisey response
x = x1;
y = y1;

% Function for returning canonical data
get_canonical_data = @(x,y)( struct('T',F(x)'*F(x), ... % Form T matrix
				    'b',F(x)'*y,    ... % Form B'y
				    'V',y'*y,       ... % Form sum of squares for y
				    'n',size(x,1)));    % The number of data points
info1 = get_canonical_data(x,y);
info = info1;

% Function for transforming canonical information into estimates
get_estimates = @(info) ( struct('ahat',info.T \ info.b, ...
				 'cov_ahat',sigma2*inv(info.T), ...
				 'hat_cov_ahat', (info.V - info.b'*inv(info.T)*info.b)/(info.n - size(info.b,1))*inv(info.T),... 
				 'hat_sigma2', (info.V - info.b'*inv(info.T)*info.b)/(info.n - size(info.b,1))));
estimates = get_estimates(info);
figure(1)
plot_estimates(info,estimates,x,y,F,a,sigma2); % See plot_estimates.m
title(sprintf('Estimates for %d points',n1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate n2 new data points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input(sprintf('Strike Enter to add %d points',n2)) 
x2 = rand(n2,1);% Generate n1 points in [0,1]
y2 = F(x2)*a + sqrt(sigma2)*randn(n2,1);% Generate noisey response
x = [x1;x2];
y = [y1;y2];

info2 = get_canonical_data(x2,y2);

% Function for combining canonical information.
combine_information = @(info1,info2)(struct('T',info1.T+info2.T,... 
                                            'b',info1.b+info2.b,    ... 
                                            'V',info1.V+info2.V,       ...  
					    'n',info1.n+info2.n));

info = combine_information(info1,info2);
estimates = get_estimates(info);

figure(2)
plot_estimates(info,estimates,x,y,F,a,sigma2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate n3 new data points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input(sprintf('Strike Enter to add %d points',n3)) 
x3 = rand(n3,1);% Generate n1 points in [0,1]
y3 = F(x3)*a + sqrt(sigma2)*randn(n3,1);% Generate noisey response
x = [x;x3]; % Although it is not necessary to keep all of the data, 
y = [y;y3]; % we do so to illustrate in the following plots

% Having defined all of the functions, the algorithm for updating information is concisely presented below
info3 = get_canonical_data(x3,y3);
info = combine_information(info,info3);
estimates = get_estimates(info);
figure(3)
plot_estimates(info,estimates,x,y,F,a,sigma2);
