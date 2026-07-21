function ANS = experfc(z)
% EXPERFC(NU,Z) calculates sqrt(Z)*exp(Z)*erf(sqrt(Z))
%
%    Matlab builtin function BESSELK(NU,Z) is identically zero for
%    Z > 700. Thus EXPBESSELK uses the asymptotic series in A&S (7.1.23)
LIMIT = 100;

ANS = zeros(size(z));

I = (z<LIMIT);
if nnz(I)
  tmp = sqrt(z(I));
  ANS(I) = tmp.*exp(z(I)).*erfc(tmp);
end

I = (z>=LIMIT);
if nnz(I)

  k = 0;
  tmp = 1/sqrt(pi);
  ANS(I) = tmp;

  while norm(tmp,inf) > 1000*eps
    k = k + 1;
    tmp = -tmp.*(2*k-1)./(2*z(I));
    ANS(I) = ANS(I) + tmp;
  end

end
