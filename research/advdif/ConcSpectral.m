function [Sig, r, th, C, H] = ConcSpectral(Pe,Nr,Nth)
% CONCSPECTRAL solves the adv-diff eqn with the spectral method.
%
%    [Sig,r,th,C,H] = ConcSpectral(Pe,Nr,Nth)
%
%    Eqn to solve:  C(r,th) = e^{eta(r,th)} * H(r,th)
%          where eta(r,th) = -log(r)/2 + (Pe/2)(2-1/r-r)(1-cos(th))
%
%    r^3 C_rr + [r^2 + Pe r(1-r^2) cos(th)] C_r 
%     + r C_thth + Pe (1+r^2) sin(th) C_th = 0
%
%    r^3 H_rr + Pe(r-r^3) H_r + r H_thth + 2Pe r sin(th) H_th
%     + (Pe r cos(th) - Pe + r/4) H = 0
%
%    See also CONCADPMESH.
%
%    @copyright
%    J. Choi, D. Margetis, T. M. Squires, and M. Z. Bazant
%    J. Fluid Mech. 2005     (cond-mat/0403740)


%%%%
%%%% Differentiation Matrices w.r.t R
%%%%

[Dr1,r] = mycheb(Nr+1,[0 1]); r = r'; % r: [0 1] (row)
Dr2 = Dr1^2;

r_ = r(2:end-1); % r_ is a (1,Nr) vector, r_: (0 1)

Dr1_endrow = Dr1(end,:); % used to get H_r(r=1)
Dr1_endcol = Dr1(2:end-1,end); % r: [0 1] --> (0 1)
Dr1 = Dr1(2:end-1,2:end-1); % Nr1 is a (Nr,Nr) matrix.

Dr2_endcol = Dr2(2:end-1,end); % r: [0 1] --> (0 1)
Dr2 = Dr2(2:end-1,2:end-1); % Nr2 is a (Nr,Nr) matrix.

%%%%
%%%% Differentiation Matrices w.r.t. TH
%%%%

[Dth1,th] = mycheb(2*Nth-2,[0 2*pi]); % th: [0 2pi]

th_ = th(1:Nth); % th_ is a (Nth,1) vector, th_: [0 2pi]

Dth2 = Dth1^2;

Dth1 = Dth1(1:Nth,1:Nth) + Dth1(1:Nth,end:-1:Nth); % th: [0 pi] --> [0 2pi]
Dth1(:,end) = Dth1(:,end)/2; 
Dth1(end,:) = 0; % It's true already because of symmetry of Dth1.
% Nth1 is a (Nth,Nth) matrix.

Dth2 = Dth2(1:Nth,1:Nth) + Dth2(1:Nth,end:-1:Nth); % th: [0 2pi] --> [0 pi]
Dth2(:,end) = Dth2(:,end)/2; 
% Nth2 is a (Nth,Nth) matrix.

%%%%
%%%% Coefficients in the P.D.E. for H(r,th)
%%%%

onesNth = ones(size(th_));
onesNr = ones(size(r_)); % not used much
eyeNth = spdiags(onesNth, 0, Nth, Nth);
eyeNr  = spdiags(onesNr', 0, Nr, Nr);
Nall = Nr*Nth;

% coef of H_rr: r^3
h_r2 = onesNth * r_.^3;
% coef of H_r: Pe(r-r^3)
h_r1 = onesNth * Pe*(r_ - r_.^3);
% coef of H_thth: r 
h_th2 = onesNth * r_;
% coef of H_th: 2 Pe r sin(th)
h_th1 = 2*Pe*sin(th_)*r_;
% coef of H: Pe r cos(th) + r/4 - Pe
h = (Pe*cos(th_) + 0.25)*r_ - Pe;

%%%%
%%%% Neumann boundary condition: H_th = 0 at th = 0
%%%% 

% looks like
% 0 H_rr + 0 H_r + 0 H_thth + 1 H_th + 0 H = 0
h_r2(1,:) = 0;
h_r1(1,:) = 0;
h_th2(1,:) = 0;
h_th1(1,:) = 1;
h(1,:) = 0;

%%%%
%%%% The linear eqn to solve:  A H = -B 
%%%%

% A is a sparse (Nr*Nth, Nr*Nth) matrix.
A = spdiags(h_r2(:),0,Nall,Nall) * kron(Dr2, eyeNth) ...
    + spdiags(h_r1(:),0,Nall,Nall) * kron(Dr1, eyeNth) ...
    + spdiags(h_th2(:),0,Nall,Nall) * kron(eyeNr, Dth2) ...
    + spdiags(h_th1(:),0,Nall,Nall) * kron(eyeNr, Dth1) ...
    + spdiags(h(:),0,Nall,Nall);


% B is a (Nr*Nth, 1) vector.
% These terms are the contribution that the boundary H(r=1,th) = 1
% make to the differenciation matrices (w.r.t. r).
B = h_r2(:) .* kron(Dr2_endcol, onesNth) ...
    + h_r1(:) .* kron(Dr1_endcol, onesNth);

%%%%
%%%% Direct Solving by '\' operator in matlab
%%%%
H_ = A\(-B);
H_ = reshape(H_,Nth,Nr);
H = [zeros(size(th_)) H_ onesNth];

Sig = H*Dr1_endrow' - 0.5; % Sig = H_r(r=1,th) - 1/2
Sig = Sig([1:end end-1:-1:1]); 

H = H([1:end end-1:-1:1],:); % th: [0 pi] --> [0 2pi]

C = H_ .* (onesNth*sqrt(1./r_)) .* exp( Pe*sin(th_/2).^2*(2-1./r_-r_) );
C = [zeros(size(th_)) C onesNth];
C = C([1:end end-1:-1:1],:); % th: [0 pi] --> [0 2pi]
