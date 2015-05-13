% FinSupFun_Demo_Sym

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
