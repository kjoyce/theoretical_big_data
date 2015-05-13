clear all, close all;

% Specify sigma^2
s = .1;

% Specify a,b, and support for r, [-T,tau]
a  = FinSupFun(.2*exp(-.1*(0:40)),0);
b  = FinSupFun(.03*(7:-1:1));
T = 7;
tau = 5;

% Construct xhat on 100 pixels
mu = FinSupFun(randn(1,100 + (a.r-a.l) + (b.r-b.l) + (T+tau)),0);
x = b .* mu;
y0 = a .* x;
y = y0 + FinSupFun(s*randn(size(y0.f)),y0.l);

% Construct P matrix
phi = b*b'; % f in the notes
p = a*phi*a' + FinSupFun(s^2); % p vector from completing the square 
q = phi*a'; % q vector from completing the square
q_delta = q.restricted_to(-T,tau); % Restrict (or zero pad) to the support [-d,d]. See FinSupFun.
r = p \ q_delta; % Solving incorported in the FinSupFun class
xhat0 = r.*y; % Construct estimate

% Calculate H
H = phi.f(-phi.l+1) + q_delta.f*r.f'; % Theoretical H = f_0 + <q,P^(-1)q>

figure(1)
hb = plot(b.l:b.r, b.f,'g','LineWidth',3);
hold on
hf = plot(phi.l:phi.r, phi.f,'m','LineWidth',3);
ha = plot(a.l:a.r, a.f,'b','LineWidth',3);
hp = plot(p.l:p.r, p.f,'r','LineWidth',3);
hold off
grid on
legend([hb hf ha hp], 'b', '\phi', 'a',  'p');

figure(2);
subplot(311)
hm = plot(mu.l:mu.r, mu.f,'g','LineWidth',3);
hold on
hx = plot(x.l:x.r, x.f,'b','LineWidth',3);
hy = plot(y.l:y.r, y.f,'r','LineWidth',3);
hold off
grid on
legend([hm hx hy], '\mu', 'x',  'y', 'Location','northwest');

subplot(312)
hold all 
plot(x.l:x.r,x.f,'LineWidth',3);
plot(xhat0.l:xhat0.r,xhat0.f,'LineWidth',3) 
hold off
title(sprintf('r PSF support [-%d,%d]',T,tau))
legend({'x',sprintf('$\\widehat{x}$, [-%d,%d]',T,tau)},'Interpreter','Latex','Location','northwest')
grid on

% Now do reconstructions for various T and tau
[T,tau] = meshgrid(1:10,1:10);
H = zeros(size(T)); % Initialize H
xhat(length(T(:))) = FinSupFun(zeros(1,100)); % Initialize struct array of reconstructions
for k = 1:length(T(:)) 
  [i,j] = ind2sub(size(T),k);
  q_delta = q.restricted_to(-T(i,j),tau(i,j)); % Restrict P and q to the support [-d,d]
  r = p \ q_delta; % Solving incorported in the FinSupFun class
  if length(r.f) < length(y.f) % Inner convolution is not symmetric. Need the smaller one to go first. 
    xhat(k) = r.*y; % Construct estimate
  else
    xhat(k) = y.*r; % Construct estimate
  end
  H(i,j) = phi.f(-phi.l+1) - q_delta.f*r.f'; % Theoretical H = f_0 + <q,P^(-1)q>
end

subplot(313)
hold all 
plot(x.l:x.r,x.f,'LineWidth',3);
plot_idx = 1:ceil(length(xhat)/4):length(xhat); % Plot 4 reconstructions
for j =  plot_idx 
  plot(xhat(j).l:xhat(j).r,xhat(j).f,'LineWidth',3) 
end
hold off
title('Reconstructions for various d')
legend(['x',arrayfun(@(idx)( sprintf('$\\widehat{x}$, [-%d,%d]',T(idx),tau(idx))),plot_idx,'UniformOutput',0)],'Interpreter','Latex','Location','northwest')
grid on

figure(3)
surf(T,tau,log(H)), title('Natural Log of Sum of Squared Error as a function of d, log H(d)'), xlabel('T'), ylabel('\tau');

%figure(4)
%idx = find(T(:) == 7 & tau(:) == 5); 
%hold all 
%plot(x.l:x.r,x.f,'LineWidth',3);
%plot(xhat(idx).l:xhat(idx).r,xhat(idx).f,'LineWidth',3) 
%hold off
%title(sprintf('r PSF support [-%d,%d]',T,tau))
%legend({'x',sprintf('$\\widehat{x}$, [-%d,%d]',T,tau)},'Interpreter','Latex')
%grid on
