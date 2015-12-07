function [beta, betatol, R2, residuals, dw] = simple_regression(data, X, significance)

% function:     simple_regression
% descrip:      Returns various statistics as a result of a simple OLS regression. 
%               The two-sided significance is 5% by default.
%
% inputs:       /data/      [n x k] doubles where n = # data points, k = # columns
%               /X/         [n x r] regressor, r = # of regressors (including scalar offset)
%
% outputs:      /beta/      coefficient estimates for the OLS
%               /betatol/   coeff stddev as a function of significance
%               /R2/        the R2 stat for the residuals by column 
%               /residuals/ the raw residuals from the OLS
%               /dw/        Durbin-Watson statistic
%

% handle args
if nargin < 3,
    significance = 0.05;
end

% extract the regression coeffs and residuals
beta = X \ data;
residuals = data - X * beta;

% calculate all residuals
%  Note: (this implementation is a memory hog, may consider col by col version later)
R2       = 1 - diag(residuals'*residuals)./diag(data'*data);
normbeta = tinv(1-0.5*significance, size(X,1)-2) * diag( sqrtm(inv(X'*X)) );
stdres   = std(residuals);
betatol  = (stdres' * normbeta')';

dw = zeros(2,size(residuals,2));
for i=1: size(residuals,2),
    [P, DW] = dwtest(residuals(:,i), X, 'approximate');
    dw(:,i) = [DW; P];
end


end  % function

