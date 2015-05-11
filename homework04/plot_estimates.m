%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              % 
% plot_estimates.m                                                             % 
%                                                                              % 
% This function plots curve fitting estimates from canonical information and   % 
% data generated in the file problem2.m                                        % 
%                                                                              % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = plot_estimates(x,n,n2,xhat,Q,y) 
  if nargin < 6
    y = NaN;
  end
  p1=plot(1:n,x,'k-', 'linewidth',2); % plot signal 
  set(gca(),'NextPlot','add'); % Hold on
  if not(isnan(y))
    plot((1-n2):(n+n2),y,'b-', 'linewidth',2); % plot measurement 
  end
  p2=plot(1:n,xhat,'r-','linewidth',2);
  p3=patch([1:n n:-1:1]', [xhat+sqrt(diag(Q));flipud(xhat-sqrt(diag(Q)))], [1,0,0], 'facealpha',.3,'linestyle',':');
  if isnan(y)
    legend([p1,p2,p3],{'Signal x','Estimate $\widehat{x}$','One Standard Deviation'},'Interpreter','Latex'); 
  else
    legend({'Signal x','Measurement y','Estimate $\widehat{x}$','One Standard Deviation'},'Interpreter','Latex'); 
  end
  title('Signal Measurement and Reconstruction')
  set(gca(),'NextPlot','replace'); % Hold off
end
