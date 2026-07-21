function ANS = expbesselk(nu,z)
% EXPBESSELK(NU,Z) calculates EXP(Z)*BESSELK(NU,K)
%
%    Matlab builtin function BESSELK(NU,Z) is identically zero for
%    Z > 700. Thus EXPBESSELK uses the asymptotic series in A&S (9.7.2)
LIMIT = 100;

ANS = zeros(size(z));

I = (z<LIMIT);
if nnz(I)
  ANS(I) = exp(z(I)).*besselk(nu,z(I));
end

I = (z>=LIMIT);
if nnz(I)
  mu = 4*nu^2;

  k = 0;
  tmp = sqrt(pi/2./z(I));
  ANS(I) = tmp;

  while norm(tmp,inf) > 1000*eps
    k = k + 1;
    tmp = tmp.*(mu-(2*k-1)^2)./(k*8*z(I));
    ANS(I) = ANS(I) + tmp;
  end

end
