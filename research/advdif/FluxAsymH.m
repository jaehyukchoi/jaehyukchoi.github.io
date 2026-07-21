function [Sig, Fn] = FluxAsymH(Pe,th,n_term);

%FLUXASYMH   Approximation of the maximum flux for Pe>>1
%
% function [Sig, Fn] = FluxAsymH(Pe,th,n_term);
%
%    Pe: scalar, th: column
%    n_term: from 0th leading order to n-th correction terms
%
%    Sig: Sig(th)
%
%    Fn = Fn( Pe(1+cos(th) ) for n = 2k,
%    Fn = Fn( Pe(1-cos(th) ) for n = 2k+1.
%
%
%    @copyright
%    J. Choi, D. Margetis, T. M. Squires, and M. Z. Bazant
%    J. Fluid Mech. 2005     (cond-mat/0403740)

if nargin<3, n_term = 100; end
if nargin<2, error('It needs at least 2 arguments'), end
if length(Pe)>1, error('Pe should be a scalar.'), end
if size(th,2)>1, error('th sould be a column vector.'), end

one_th = ones(size(th));
n_th = length(th);

if n_term == 0

  Fn = one_th;
  Sig = 2*sqrt(Pe/pi)*abs(sin(th/2));
  
else

  cos_th = cos(th);
  u = Pe*(1-[-cos_th; cos_th]); % even  odd
  [u_,tmp,ind] = unique(u);
  one_u = ones(size(u_));
  n_u = length(u_);
  
  [t, I] = myclencurt(200, [0 inf]); %z: integration variable [0 pi/2]
  t2 = t.^2;
  one_t = ones(size(t));
  
  [t2x,t2y] = meshgrid(t2);
  R_I = repmat(t2'.*I, size(t2))./(2 + t2x + t2y); % integration along the row
  Q = one_u * [(1/pi)*exp(-2*Pe*t2)./sqrt(2+t2)];
  
  Fn = zeros(n_u, n_term);

  % k = 1
  tmp = sqrt(1/pi) * one_u*[exp(-2*Pe*t2)./sqrt(2+t2)].* experfc(u_*(2+t2));
  fn = Q - tmp;
  Fn(:,1) = (1/pi)*expbesselk(0,2*Pe) - tmp*I;
  
  for k = 2:n_term
    fn = Q .* ( fn * R_I );  % Q int R(t2,t1) fn(t1) dt1 
    Fn(:,k) = fn*I;  % int fn(t) dt
    if nargin <3 && max(Fn(:,k)) < 1000*eps, 
      n_term = k; break;
    end
  end

  Fn = Fn(ind,1:n_term);
  % odd
  Fn_odd = Fn(n_th+1:end,1:2:n_term);
  Fn = Fn(1:n_th,:);
  Fn(:,1:2:n_term) = Fn_odd;
  Fn = [one_th Fn];

  % Construct prefactor:
  %
  % PreF = e^{-2n Pe} |sin(th/2)|  for  n = 2k
  %        e^{-2(n-cos(th)) Pe} |cos(th/2)|  for n = 2k+1
  %
  
  PreF = zeros(n_th,n_term+1);
  
  n_even = [0:2:n_term];
  PreF(:,n_even+1) = abs(sin(th/2)) * exp(-2*Pe*n_even);
  n_odd = [1:2:n_term];
  PreF(:,n_odd+1) =  [abs(cos(th/2)).*exp(2*Pe*cos(th))] * exp(-2*Pe*n_odd);

  Sig = 2*sqrt(Pe/pi) .* sum(PreF.*Fn,2);

end
