function [x,w] = myclencurt(N,xlim)
% MYCLENCURT computs weights for Clenshaw-Curtis quadrature
%
%    [x,w] = myclencurt(N, xlim)
%
%    x: grid points (row), w: weights (column),
%    xlim = [xi xf] or [0 1] by default.
%           [0 inf] or [-inf inf] is also acceptable.
%
%    Modified from the original m-file "clencurt.m" in
%    "Spectral Methods in Matlab" by L. N. Trefethen (Siam, 2000)
%    http://web.comlab.ox.ac.uk/oucl/work/nick.trefethen/

if mod(N,2)~=0, error('N should be an even integer.'), end
if nargin < 2, xlim = [0 1]; end

%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%
% orignal m-file: clencurt.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = pi*(0:N)'/N; x = cos(theta);
w = zeros(1,N+1); ii = 2:N; v = ones(N-1,1);
if mod(N,2)==0 
  w(1) = 1/(N^2-1); w(N+1) = w(1);
  for k=1:N/2-1, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
  v = v - cos(N*theta(ii))/(N^2-1);;
else
  w(1) = 1/N^2; w(N+1) = w(1);
  for k=1:(N-1)/2, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
end
w(ii) = 2*v/N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = w'; % column
x = x'; % row

if xlim(2) == inf
  
  if xlim(1) == -inf  % case: xlim = [0 inf]
    x = tan(pi/2*x(end-1:-1:2));
    w = (pi/2) * (1+x'.^2).*w(end-1:-1:2);
  
    
  elseif xlim(1) == 0 % case: xlim = [-inf inf]
    x = tan(pi/2*x(N/2+1:-1:2));
    w = (pi/2) * (1+x'.^2).*w(N/2+1:-1:2);
    w(2:end) = 2*w(2:end);
  
  else
    error('Make sure that xlim = [0 inf] or [-inf inf]')
  end

  
else % case = [finite finite]
  
  x = xlim(1)*(1/2+x/2) + xlim(2)*(1/2-x/2);
  w = w*abs(diff(xlim))/2;

end
