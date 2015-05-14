clear all, close all;

% Specify sigma^2
s = .1;

% Specify a,b, and r width
a = FinSupFun(.1+0*(1:15));
b = FinSupFun(.03*(7:-1:1));
d = 120;

% Construct xhat on 100 pixels
mu = FinSupFun(randn(1,100 + (a.r-a.l) + (b.r-b.l) + 2*d),0);
x = b .* mu;
y0 = a .* x;
y = y0 + FinSupFun(s*randn(size(y0.f)),y0.l);

% Construct P matrix
phi = b*b'; % f in the notes
p = a*phi*a' + FinSupFun(s^2); % p vector from completing the square 
q = phi*a'; % q vector from completing the square
q_delta = q.restricted_to(-d,d); % Restrict (or zero pad) to the support [-d,d]. See FinSupFun.
r = p \ q_delta; % Solving incorported in the FinSupFun class
xhat0 = r.*y; % Construct estimate
del = r*a; % Estimate of delta function

% Calculate H
H = phi.f(-phi.l+1) + q_delta.f*r.f'; % Theoretical H = f_0 + <q,P^(-1)q>

figure(1)
subplot(121)
hb = plot(b.l:b.r, b.f,'g','LineWidth',3);
hold on
hf = plot(phi.l:phi.r, phi.f,'m','LineWidth',3);
ha = plot(a.l:a.r, a.f,'b','LineWidth',3);
hp = plot(p.l:p.r, p.f,'r','LineWidth',3);
hr = plot(r.l:r.r, r.f,'c','LineWidth',3);
hd = plot(del.l:del.r, del.f,'k','LineWidth',3);
hold off
grid on
legend([hb hf ha hp hr hd], {'b', '$\phi$', 'a',  'p', 'r', '$r * a$'}, 'Interpreter','Latex');
title(sprintf('Covariance PSFs for Deconvolution d = %d',d));

figure(2);
subplot(311)
hm = plot(mu.l:mu.r, mu.f,'g','LineWidth',3);
hold on
hx = plot(x.l:x.r, x.f,'b','LineWidth',3);
hy = plot(y.l:y.r, y.f,'r','LineWidth',3);
hold off
grid on
legend([hm hx hy], '\mu', 'x',  'y');

subplot(313)
hold all 
plot(x.l:x.r,x.f,'LineWidth',3);
plot(xhat0.l:xhat0.r,xhat0.f,'LineWidth',3) 
patch([xhat0.l:xhat0.r xhat0.r:-1:xhat0.l],[xhat0.f+sqrt(H) fliplr(xhat0.f-sqrt(H))], [0 1 0], 'facealpha', .3);
hold off
title(sprintf('r PSF radius d=%d',120))
legend({'x','$\widehat{x}$, d=120'},'Interpreter','Latex')
grid on

% Now do reconstructions for d from 5 to 130
d = 5:130;
H = zeros(size(d)); % Initialize H
xhat(length(d)) = FinSupFun(zeros(1,100)); % Initialize struct array of reconstructions
r(length(d)) = FinSupFun(zeros(1,100)); % Initialize struct array of reconstructions
for i = 1:length(d) 
  q_delta = q.restricted_to(-d(i),d(i)); % Restrict P and q to the support [-d,d]
  r(i) = p \ q_delta; % Solving incorported in the FinSupFun class
  if length(r(i).f) < length(y.f) % Inner convolution is not symmetric. Need the smaller one to go first. 
    xhat(i) = r(i).*y; % Construct estimate
  else
    xhat(i) = y.*r(i); % Construct estimate
  end
  H(i) = phi.f(-phi.l+1) - q_delta.f*r(i).f'; % Theoretical H = f_0 + <q,P^(-1)q>
end 

subplot(312)
hold all 
plot(x.l:x.r,x.f,'LineWidth',3);
plot_idx = 1:ceil(length(xhat)/4):length(xhat); % Plot 4 reconstructions
for j =  plot_idx 
  plot(xhat(j).l:xhat(j).r,xhat(j).f,'LineWidth',3) 
end
hold off
title('Reconstructions for various d')
legend(['x',arrayfun(@(d)( sprintf('$\\widehat{x}$, d=%d',d)),d(plot_idx),'UniformOutput',0)],'Interpreter','Latex')
grid on

figure(1)
subplot(122)
title('PSFs for a*f')
hold all
for j = plot_idx
  del = r(j) * a;
  plot(del.l:del.r,del.f, 'linewidth',3); 
end
hold off
grid on
legend([arrayfun(@(d)( sprintf('$r_d * a$, d=%d',d)),d(plot_idx),'UniformOutput',0)],'Interpreter','Latex')

figure(3)
plot(d,H,'ko'), title('Sum of Squared Error as a function of d, H(d)'), xlabel('d'), ylabel('H(d)');

