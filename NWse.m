function VNW = nw(x,m,kernel,err)
% See Pesaran (2015) for formulas and notations
% Input:        - x: T*n matrix of independent variables
%               - err: errors of regression
%               - m: bandwidth, i.e. number of lags
%               - kernel: 1 = uniform, 2 = Bartlett (default)
% Output:       - Newey-West Variance matrix
% If no errors are provided, computation is different
% Written by Lars Spreng

no_err = 0;
if nargin < 4
    no_err = 1;
    err = ones(100,1);
    display('No regression errors were provided; proceeding without')
end
if nargin < 3
    kernel = 2;
end
T = size(x,1);
Q = x'*x./T;
Omega_0 = ((x.*err)'*(x.*err));
for jj=1:m
    
    if kernel == 1
        w = 1;
    elseif kernel == 2
        w =  1 - (jj/(m+1));
    end
    
    Omega = (x(jj+1:end,:).*err(jj+1:end))'*(x(1:end-jj,:).*err(1:end-jj));
    Omega_0 = Omega_0 + w*(Omega + Omega');
end
S = Omega_0;

if no_err == 0
    VNW = 1/T*inv(Q)*S*inv(Q);
elseif no_err == 1
    VNW = 1/T*S;
end
