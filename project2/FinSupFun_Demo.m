% FinSupFun_Demo
clear all; close all;
s = .1;

mu = FinSupFun(randn(1,200),0);
b  = FinSupFun(.03*(7:-1:1));
x  = b .* mu;
a(1)  = FinSupFun(.1+0*(1:15));
a(2)  = FinSupFun(.2*exp(-.1*(0:40)),0);

y0 = a(1) .* x;
y(1)  = y0 + FinSupFun(s*randn(size(y0.f)), y0.l);
phi = b*b';
p(1) = a(1)*phi*a(1)' + FinSupFun(s^2);

y0 = a(2) .* x;
y(2)  = y0 + FinSupFun(s*randn(size(y0.f)), y0.l);
phi = b*b';
p(2) = a(2)*phi*a(2)' + FinSupFun(s^2);

figure(1);
subplot(211)
hm = plot(mu.l:mu.r, mu.f,'g','LineWidth',3);
hold on
hx = plot(x.l:x.r, x.f,'b','LineWidth',3);
hy = plot(y(1).l:y(1).r, y(1).f,'r','LineWidth',3);
hold off
grid on
legend([hm hx hy], '\mu', 'x',  'y');
title('Symmetric blur by a')
subplot(212)
hm = plot(mu.l:mu.r, mu.f,'g','LineWidth',3);
hold on
hx = plot(x.l:x.r, x.f,'b','LineWidth',3);
hy = plot(y(2).l:y(2).r, y(2).f,'r','LineWidth',3);
hold off
grid on
legend([hm hx hy], '\mu', 'x',  'y');
title('Asymmetric blur by a')

figure(2);
subplot(121)
hb = plot(b.l:b.r, b.f,'g','LineWidth',3);
hold on
hf = plot(phi.l:phi.r, phi.f,'m','LineWidth',3);
ha = plot(a(1).l:a(1).r, a(1).f,'b','LineWidth',3);
hold off
grid on
legend([hb hf ha ], 'b', '\phi', 'a');
title('Symmetric blur by a')
subplot(122)
hb = plot(b.l:b.r, b.f,'g','LineWidth',3);
hold on
hf = plot(phi.l:phi.r, phi.f,'m','LineWidth',3);
ha = plot(a(2).l:a(2).r, a(2).f,'b','LineWidth',3);
hold off
grid on
legend([hb hf ha ], 'b', '\phi', 'a');
title('Asymmetric blur by a')
