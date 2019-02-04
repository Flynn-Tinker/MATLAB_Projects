function [rho,rL,rU] = kb_bootstrap_rho(x,y,varargin)
% kb_bootstrap_rho(x,y,alpha,N) does a very basic bootstrap calculation to 
%   determine confidence intervals on the pearson correlation.  Optional 
%   arguments are the desired confidence limit (defaults to 0.05), and the 
%   number of bootstrap replications (defaults to 2000).  CIs are computed 
%   using the (simple) interval method; Efron's accelerator is likely better 
%   but more complicated.
%
%       Example:
%           [rho,rL,rU] = kb_bootstrap_rho(x,y,'p',0.05)
%           [rho,rL,rU] = kb_bootstrap_rho(x,y,'N',2000)
%           [rho,rL,rU] = kb_bootstrap_rho(x,y,'p',0.05,'N',2000)
%       
%       Input:
%           x,y : vectors for which the correlation is desired. (If x or
%                 y has significant autocorrelation you need a more
%                 complicated version of the bootstrap).
%           p   : (optional) significance limit; defaults to 0.05 (so 
%                 the returned (rL,rU) will cover 95% of the bootstrap
%                 distribution
%           N   : (optional) number of bootstrap replications; defaults
%                 to 2000.
%
%       Output:
%           rho  : the actual observed correlation for the real x,y
%           rL   : lower limit of the CI from the bootstrap
%           rU   : upper limit of the CI from the bootstrap
%

% setup default values
pars = inputParser;
pars.addRequired('x',@isfloat);
pars.addRequired('y',@isfloat);
pars.addParamValue('p',0.05,@isfloat);
pars.addParamValue('N',2000,@isnumeric);

% parse
pars.parse(x,y,varargin{:});

% assign arguments
x = pars.Results.x(:);
y = pars.Results.y(:);
p = pars.Results.p;
N = pars.Results.N;

% observed correlation
rho = corr(x,y);

% bootstrap
nSamp = size(x,1);

rhobs = zeros(N,1);
for iB = 1:N
    % draw with replacement - works with or without the stats toolbox
    try
        randx = unidrnd(nSamp,nSamp,1);
    catch
        randx = randint(nSamp,1,nSamp);
    end
    % compute the correlation
    rhobs(iB) = corr(x(randx),y(randx));
end

% sort the values to obtain the correlation coefficient
y = sort(rhobs);
rL = y(floor(N*p/2));
rU = y(floor(N*(1-p/2)));

% hist(rhobs,25)

return