function [r,r_x1,r_x2] = AdpmeshR(Pe, x)
% ADPMESHR generates adaptive mesh for radial direction.
% 
%     [r,r_x1,r_x2] = AdpmeshR(Pe, x)
%
%     ADPMESHR and ADPMESHTH are used in FLUXADPMESH

%%%%%% Another Mesh Using Essential Singularity

LOW = 0.1; HIGH = 10;

if Pe < LOW | HIGH < Pe

  HighPe = 0;
  if HIGH < Pe, Pe = 1/sqrt(Pe); x = 1-x; HighPe = 1; end

  a = (1-Pe);
  r = x - (a/pi)*sin(pi*x);
  r_x1 = 1 - a*cos(pi*x);
  r_x2 = (a*pi)*sin(pi*x);
    
  if HighPe, r = 1-r; r_x2 = -r_x2; end

else
  r = x;
  r_x1 = ones(size(x));
  r_x2 = zeros(size(x));
end

%%%% ANOTHER ADAPTIVE MESH
%  a = 8;   % increase 'a' if you want more stretch at the b-layer
%  xo = Pe; for n=1:5, xo = a / ( a - log(Pe/(1+a*(1-xo)/xo^2)) ); end
%  ro = exp(a*(1-1/xo));
%  x = (1-xo).*x + xo;
%  f = exp( a*(1-1./x) );
%  r = (f-ro)/(1-ro);
%  r_x1 = (1-xo)/(1-ro) * (a./x.^2).*f;
%  r_x2 = (1-xo)^2/(1-ro) * (a^2./x.^4 - 2*a./x.^3).*f;
   
