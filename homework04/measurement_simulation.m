%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
% measurement_simulation.m                                                     %
%                                                                              %
% The following codes simulates signal estimation from various measurement     %
% models following the methods outlined in Lectures 08a-08c of Math 491 Spring %
% 2015.                                                                        %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem 1
clear all; %close all;
n = 100;	  % Number of sample points of x
n1 = 7;		  % Kernel radius of F^(1/2) 
n2 = 9;           % Kernel radius for Blurring operator A
s = .1;	          % Standard Deviation of measurment noise

% Definitions for various kernels
ker.rectangular = @(d)( ones(1,2*d+1) ); 
ker.triangular  = @(d)( [1:d+1 d:-1:1] );
ker.parabolic   = @(d)( (d+1)^2-(-d:d).^2 );
ker.gauss       = @(d)( exp(-(-d:d).^2/(2*(d/3)^2)) );

% Kernel for the prior F^(1/2)
b0 = ker.triangular(n1); b0 = b0(((size(b0,2)+1)/2):end); b = [b0 zeros(1,n-size(b0,2))];
B = toeplitz(b); % prior F^(1/2)
F = B'*B;        % prior variance

x = B*randn(size(F,1),1);  % Generate x 

a = ker.rectangular(n2)/sum(ker.rectangular(n2)); % Kernel for Blurring operator A
A = convmtx(a',n);  % Create convolution matrix

y = A*x + s*randn(size(A,1),1); % Generate noisey data
figure(1), plot(1:n,x,'k-', (1-n2):(n+n2),y,'b-', 'linewidth',2); % plot signal and measurement
title('Signal and Measurement');
grid on;

% Problem 2
Q = inv(1/s*A'*A + inv(F)); % Estimator variance
xhat = Q*(1/s*A'*y); % Estimator
figure(2), plot_estimates(x,n,n2,xhat,Q,y); grid on;
title('Signal, Measurement, and Reconstruction with Prior Information');

% Problem 3 (a)
get_canonical_data = @(y,A,S)( struct('T',A'*inv(S)*A, ... % Form T matrix
				      'b',A'*inv(S)*y  ... % Form B'y
				    )); 
info1 = get_canonical_data(y,A,s*eye(n+2*n2));
get_estimates = @(info)( struct('Q',inv(info.T), ...
                                'xhat',(info.T)\info.b ...
				 ));
estimates1 = get_estimates(info1);
figure(3), plot_estimates(x,n,n2,estimates1.xhat,estimates1.Q,y); grid on;
title('Signal, Measurement, and Reconstruction with No Prior Information');

% Problem 3 (b)
combine_information = @(info1,info2)(struct('T',info1.T+info2.T,...
					    'b',info1.b+info2.b));

prior_to_canonical = @(mu,F)( get_canonical_data(mu,eye(n),F) );
prior_info = prior_to_canonical(zeros(n,1),F);
info2 = combine_information(info1,prior_info);
estimates2 = get_estimates(info2);
figure(4), plot_estimates(x,n,n2,estimates2.xhat,estimates2.Q,y); grid on;
title('Signal, Measurement, and Reconstruction with Prior Information');

% Problem 3 (c)
figure(5);
plot(1:n,x,'k-', 'linewidth',2); % plot signal 
hold all;
names = fieldnames(ker);
info = get_canonical_data(zeros(n,1),zeros(n),eye(n)); % No information
for j = 1:length(names);
  a = ker.(names{j})(n2)/sum(ker.(names{j})(n2));
  a = a/sum(a);
  A = convmtx(a',n);  % Create convolution matrix 
  y = A*x + s*randn(size(A,1),1); % Generate noisey data
  new_info = get_canonical_data(y,A,s*eye(n+2*n2));
  info = combine_information(info,new_info);
  plot((1-n2):(n+n2),y);
end 
hold off;
grid on;
legend(['true signal x';fieldnames(ker)]);
title('Measurements from 4 different blurring kernels');

estimates = get_estimates(info);
figure(6), plot_estimates(x,n,n2,estimates.xhat,estimates.Q);
title('Reconstruction from 4 different blurring kernels');

% Multiple measurements
figure(7);
plot(1:n,x,'k-', 'linewidth',2); % plot signal 
hold all;
names = fieldnames(ker);
info = get_canonical_data(zeros(n,1),zeros(n),eye(n)); % No information
nmeas = 10;
for jj = 1:(nmeas*length(names));
  j = mod(jj,length(names)) + 1;
  a = ker.(names{j})(n2)/sum(ker.(names{j})(n2));
  a = a/sum(a);
  A = convmtx(a',n);  % Create convolution matrix 
  y = A*x + s*randn(size(A,1),1); % Generate noisey data
  new_info = get_canonical_data(y,A,s*eye(n+2*n2));
  info = combine_information(info,new_info);
  h = plot((1-n2):(n+n2),y);
end 
hold off;
grid on;
legend(['true signal x';fieldnames(ker)]);
title('10 measurements each from 4 blurring kernels $(n=40)$','Interpreter','Latex');

estimates = get_estimates(info);
figure(8), plot_estimates(x,n,n2,estimates.xhat,estimates.Q);
title('Reconstruction from 10 measurements each from 4 blurring kernels $(n=40)$','Interpreter','Latex')

% Problem 3 (d)
names = fieldnames(ker);
info = prior_to_canonical(zeros(n,1),F);
for j = 1:length(names);
  a = ker.(names{j})(n2)/sum(ker.(names{j})(n2));
  a = a/sum(a);
  A = convmtx(a',n);  % Create convolution matrix 
  y = A*x + s*randn(size(A,1),1); % Generate noisey data
  new_info = get_canonical_data(y,A,s*eye(n+2*n2));
  info = combine_information(info,new_info);
end 
estimates = get_estimates(info);
figure(9), plot_estimates(x,n,n2,estimates.xhat,estimates.Q);
title('Reconstruction from 4 different blurring kernels with Prior Information');
