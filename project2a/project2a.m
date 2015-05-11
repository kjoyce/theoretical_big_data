clear all, close all;

% Specify sigma^2
s = .1;

% Specify a,b, and r width
a = FinSupFun(.1+0*(1:15));
b = FinSupFun(.03*(7:-1:1));
d = 30;

% Always construct xhat on 100 pixels
mu = FinSupFun(randn(1,100 + (a.r-a.l) + (b.r-b.l) + 2*d),0);
x = b .* mu;
y0 = a .* x;
y = y0 + FinSupFun(s*randn(size(y0.f)),y0.l);

% Construct P matrix
phi = b*b';
p = a*phi*a' + FinSupFun(s^2);
%toeplitz_row = [p.f( (p.r+1):end ), zeros(1, mu.r-(p.r))];% This needs to be from the center of p
%P = toeplitz(toeplitz_row); 
q = phi*a';
q = q.restricted_to(-d,d);

% Restrict P and q to the support [-d,d]
%rr = P(1:(2*d+1),1:(2*d+1)) \ q.f((q.r+1-d):(q.r+1+d))';
%r = FinSupFun(rr',-d);
r = p \ q;
xhat = r.*y;

% Calculate H
%resid = xhat + FinSupFun(-x.f,x.l);
%H = norm(resid.f(d:(end-d))) 
%phi.f(-phi.l+1) + q_delta*rr % Theoretical H = f_0 + <q,P^(-1)q>

figure(1);
subplot(211)
hm = plot(mu.l:mu.r, mu.f,'g','LineWidth',3);
hold on
hx = plot(x.l:x.r, x.f,'b','LineWidth',3);
hy = plot(y.l:y.r, y.f,'r','LineWidth',3);
hold off
grid on
legend([hm hx hy], '\mu', 'x',  'y');

subplot(212)
plot(x.l:x.r,x.f,'b-','LineWidth',3) 
hold on 
plot(xhat.l:xhat.r,xhat.f,'k-','LineWidth',3) 
hold off
title(sprintf('r PSF radius d=%d',d))
legend({'x','$\widehat{x}$'},'Interpreter','Latex')
grid on

