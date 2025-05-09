function output = logdensity2loglik(logdensity,x,del,param)
% Inputs:
%
% logdensity: is the function handle of the transition density with the
% structure such as logdensity(x,x0,del,param)
%
% x: the price series nxk, n is the number of observations and k is the
% dimension of the multivariate series
%
% del: sampling interval
% param: parameter vector

n = length(x) - 1;
output = 0;
for i=1:n
    output = output + logdensity(x(i+1,:),x(i,:),del,param); 
end
