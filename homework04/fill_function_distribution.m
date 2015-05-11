%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              % 
% fill_funciton_distribution.m                                                 % 
%                                                                              % 
%  f = a matrix whose columns are functions to be represented                  % 
%  quantiles = a vector of quantiles to draw for each fill patch               % 
%                                                                              % 
%  Author: Kevin Joyce                                                         % 
%  Date: Feb 2015                                                              % 
%                                                                              % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p] = fill_funciton_distribution(f,quantiles,x,ax)
[n,m] = size(f);
if nargin < 3
  x = (1:n)';
end
if nargin < 4
  ax = gca();
end
hold(ax,'on');
for middle_percent = quantiles
  f_sorted = sort(f,2);
  top_idx = ceil(m*(1-(1-middle_percent)/2));
  bot_idx = floor(m*(1-middle_percent)/2);
  if bot_idx < 1
    bot_idx = 1;
 end
  p = patch([x;flipud(x)], [f_sorted(:,top_idx); flipud(f_sorted(:,bot_idx))], [1,0,0],'facealpha',.3,'edgecolor','none');
end
hold(ax,'off');
