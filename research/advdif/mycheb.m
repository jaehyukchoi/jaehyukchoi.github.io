function [D,x] = mycheb(N,xlim)
% MYCHEB compute Chebshev differentiation matrix.
%
%    [D,x] = mycheb(N,xlim)
%
%    D: differentiation matrix, x:grid points,
%    xlim = [xi xf] or [0 1] by default
%
%    Modified from the original m-file "cheb.m" in
%    "Spectral Methods in Matlab" by L. N. Trefethen (Siam, 2000)
%    http://web.comlab.ox.ac.uk/oucl/work/nick.trefethen/

if nargin < 2, xlim = [0 1]; end

%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%
% orignal m-file: chen.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N==0, D=0; x=1; return, end

x = cos(pi*(0:N)/N)'; 
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';                  
D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
D  = D - diag(sum(D'));                 % diagonal entries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = xlim(1)*(1+x)/2 + xlim(2)*(1-x)/2;
D = D*2/(xlim(1)-xlim(2));
