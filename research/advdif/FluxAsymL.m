function Sig = FluxAsymL(Pe,th);
%FLUXASYM    Approximation of the maximum flux for Pe<<1
%
% function Sig = FluxAsymL(Pe,th,High);
%
%    For Pe<<1, MAXFLUX(Pe) = 1./(-gamma-log(Pe/4)), where gamma = -psi(1) is
%    Euler's constant
%    For Pe>>1, MAXFLUX(Pe) = 2*sqrt(Pe/pi)
%
%
%    @copyright    
%    J. Choi, D. Margetis, T. M. Squires, and M. Z. Bazant
%    J. Fluid Mech. 2005     (cond-mat/0403740)

if nargin<2, th = pi; end
if length(th)>1 && length(Pe)>1, error('Only one of th and Pe can be a vector.'), end

if length(Pe)>1, error('Pe should be a scalar.'), end
if size(th,2)>1, error('th sould be a column vector.'), end

n_th = length(th);
one_th = ones(n_th,1);

N = 200; % You may increase this number for higher accuracy 
[t,I] = myclencurt(N,[0 Pe]); % t: integration variable [0 Pe]
F = zeros(length(th),length(I));
F(:,1) = 0.5;
t_ = t(2:end);
F(:,2:end) = exp(cos(th)*t_) .* [one_th*(besseli(1,t_)./t_)];

Sig = (besseli(0,Pe)/besselk(0,Pe/2))*exp(Pe*cos(th)) - Pe*(cos(th) + F*I);

%Sig = (besseli(0,Pe)/(psi(1)+log(4/Pe)))*exp(Pe*cos(th)) - Pe*(cos(th) + F*I);
