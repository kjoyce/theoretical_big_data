% FinSupFun_Demo_Sym
clear all, close all;
s = .1;

mu = FinSupFun(randn(1,200),0);
b  = FinSupFun(.03*(7:-1:1));
x  = b .* mu;
a  = FinSupFun(.1+0*(1:15));
y0 = a .* x;
y  = y0 + FinSupFun(s*randn(size(y0.f)), y0.l);
phi = b*b';
p = a*phi*a' + FinSupFun(s^2);

figure(1);
hm = plot(mu.l:mu.r, mu.f,'g','LineWidth',3);
hold on
hx = plot(x.l:x.r, x.f,'b','LineWidth',3);
hy = plot(y.l:y.r, y.f,'r','LineWidth',3);
hold off
grid on
legend([hm hx hy], '\mu', 'x',  'y');

figure(2);
hb = plot(b.l:b.r, b.f,'g','LineWidth',3);
hold on
hf = plot(phi.l:phi.r, phi.f,'m','LineWidth',3);
ha = plot(a.l:a.r, a.f,'b','LineWidth',3);
hp = plot(p.l:p.r, p.f,'r','LineWidth',3);
hold off
grid on
legend([hb hf ha hp], 'b', '\phi', 'a',  'p');

d = 15; % Change here
P = toeplitz(p.f( (p.r+1):end ));
q = phi*a';
rr = P(1:(2*d+1),1:(2*d+1)) \ q.f((q.r+1-d):(q.r+1+d))'; 
r = FinSupFun(rr,-(length(rr)-1)/2);
del = r*a;
xhat = r*y;
%xhat.l = 2*d + xhat.l; % Not sure why I have to do this
%xhat.r = 2*d + xhat.r; % Not sure why I have to do this
figure(1)
hold on
idx = 2*d+(xhat.l:xhat.r);
iidx = ((x.l-5)+idx(1)):((x.l-5)+idx(end));
plot(idx,xhat.f,'LineWidth',3) 

d = 3:20; 
H = zeros(1,length(d));
for i = 1:length(d)
  P = toeplitz(p.f( (p.r+1):end ));
  q = phi*a';
  rr = P(1:(2*d(i)+1),1:(2*d(i)+1)) \ q.f((q.r+1-d(i)):(q.r+1+d(i)))'; 
  r = FinSupFun(rr,-(length(rr)-1)/2);
  xhat = r.*y;
  idx = 2*d(i)+(xhat.l:xhat.r);
  iidx = ((x.l-5)+idx(1)):((x.l-5)+idx(end));
  H(i) = norm(xhat.f-x.f(iidx));
end

%q = phi * a';
%q = FinSupFun(q.f(q.l+d
%
%
%P = toeplitz(toeplitz_row);
%diff = (size(q.f,2) - size(P,2))/2;
%q_delta = q.f((diff+1):(end-diff));
%r_k = P \ q_delta';
%r = FinSupFun(r_k


