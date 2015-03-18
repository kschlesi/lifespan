function [theP,issig] = ztop(zscore,is1tailed,alpha)
% ZTOP returns the equivalent p-value of an input z-score,
% assuming a normal z-score distribution of mean 0 and std 1.
% INPUTS: zscore - a matrix of any size containing z-scores from
%                   a normal distribution with mean 0 and std 1.
%         is1tailed - a boolean value indicating whether to use a 
%                     one-tailed test (default = 0, two-tailed)
%         alpha - a value in [0,1] below which the p-value must fall to
%                 reject the null hypothesis (default = 0.05)
% OUTPUT: theP - a matrix of the same size as zscore, containing p-values
%         issig - a matrix of the same size as zscore, containing 1 if the
%                 null hypothesis is rejected and 0 if it is kept

if nargin<3
    alpha = 0.05;
end
if nargin<2
    is1tailed = 0;
end

tailfac = -1*boolean(is1tailed)+2;
theP = tailfac*(1-normcdf(abs(zscore),0,1));
issig = (theP<alpha);

end