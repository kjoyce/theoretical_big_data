% Homework #4 Tempate

M = 100;  % Dimension of Unknown Signal x
s = .1;   % Standard deviation for random noise (Var = s^2)

g = randn(M,1);

figure(1);  plot(g,'linewidth',3);
title('Uncorrelated Random Signal \mu to generate x','FontSize',14);

b0 = 7:-1:1; b = [b0 zeros(1,M-size(b0,2))];  %
B = toeplitz(b);  % Symmetric Toeplitz matrix with "triangular" profile

F = B*B'; % Variance matrix for x

x = B*g;  % Random signal with variance matrix F

figure(2);  plot(x,'linewidth',3);
title('Input Signal    x','FontSize',14);

% Profile of he Point Spread Function a
d=9; al=-d; ar=d; a = ones(1,ar-al+1); % Rectangular Profile
% d=9; al=-d; ar=d; a = [1:d+1 d:-1:1]; % Triangular Profile
% d=12; al=-d; ar=d; a = (d+1)^2-(-d:d).^2; % Parabolic Profile
% d=10; al=-d; ar=d; a = exp(-(-d:d).^2/(2*(d/3)^2)); % Gaussian Profile
% d=50; al=0; ar=d; a = exp(-(0:d)/(d/4)); % Exponential (asymmetric) Profile

a = a/sum(a);  % 'Normalized' PSF: sum a_i = 1

figure(3);  plot(al-1:ar+1,[0 a 0],'linewidth',3);
title('Point Spread Function   a','FontSize',14);

N = M + ar-al;  % Dimension of observed signal y

A = zeros(N,M);  % Matrix A determined by PSF a
for j=1:M
    A(j+(0:ar-al),j) = a;
end

y = A*x + s*randn(N,1);  % Measurement

figure(4);  plot(y,'linewidth',3);
title('Observation     y = Ax + \nu','FontSize',14);

figure(5);
plot([1+al M+ar],[0 0],'k','linewidth',2);
hold on
hx=plot(1:M,x,'linewidth',3);
hy=plot(1+al:M+ar,y,'r','linewidth',3);
title('Signal x and Observation y','FontSize',14);
legend([hx hy],'x','y');
grid on
ax=axis;  axis([1+al M+ar ax(3:4)]);

Q = inv(1/s*A'*A + inv(F));
xhat = Q*(1/s*A'*y);
plot(1:M,xhat,'k-',1:M,xhat+sqrt(diag(Q)),'k:',1:M,xhat-sqrt(diag(Q)),'k:','linewidth',2);
hold off
