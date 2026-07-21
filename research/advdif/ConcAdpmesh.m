function [Sig, r, th, C, H, rhat, thhat] = ConcAdpmesh(Pe,Nr,Nth)
% CONCADPMESH solves the adv-diff eqn with spectral method and adapeive mesh.
%
%    See also CONCSPECTRAL
%
%    @copyright
%    J. Choi, D. Margetis, T. M. Squires, and M. Z. Bazant
%    J. Fluid Mech. 2005     (cond-mat/0403740)

%%%%
%%%% Differentiation Matrices w.r.t RHAT
%%%%

[Drhat1,rhat] = mycheb(Nr+1,[0 1]); rhat = rhat'; % r: [0 1] (row)
Drhat2 = Drhat1^2;

rhat_ = rhat(2:end-1); % rhat_ is a (1,Nr) vector, rhat_: (0 1)

Drhat1_endrow = Drhat1(end,:); % used to get H_rhat(rhat=1)
Drhat1_endcol = Drhat1(2:end-1,end); % rhat: [0 1] --> (0 1)
Drhat1 = Drhat1(2:end-1,2:end-1); % Nrhat1 is a (Nr,Nr) matrix.

Drhat2_endcol = Drhat2(2:end-1,end); % rhat: [0 1] --> (0 1)
Drhat2 = Drhat2(2:end-1,2:end-1); % Nrhat2 is a (Nr,Nr) matrix.

%%%%
%%%% Differentiation Matrices w.r.t. THHAT
%%%%

[Dthhat1,thhat] = mycheb(2*Nth-2,[0 2*pi]); % thhat: [0 2pi]

thhat_ = thhat(1:Nth); % thhat_ is a (Nth,1) vector, thhat_: [0 2pi]

Dthhat2 = Dthhat1^2;

Dthhat1 = Dthhat1(1:Nth,1:Nth) + Dthhat1(1:Nth,end:-1:Nth); % thhat: [0 pi] --> [0 2pi]
Dthhat1(:,end) = Dthhat1(:,end)/2; 
Dthhat1(end,:) = 0; % It's true already because of symmetry of Dthhat1.
% Nthhat1 is a (Nth,Nth) matrix.

Dthhat2 = Dthhat2(1:Nth,1:Nth) + Dthhat2(1:Nth,end:-1:Nth); % thhat: [0 2pi] --> [0 pi]
Dthhat2(:,end) = Dthhat2(:,end)/2; 
% Nthhat2 is a (Nth,Nth) matrix.

%%%%
%%%% (R, TH) <- (RHAT, THHAT)
%%%%
[r_, r_rhat1, r_rhat2] = AdpmeshR(Pe, rhat_);
r = [0 r_ 1];
[th_, th_thhat1, th_thhat2] = AdpmeshTH(Pe, thhat_);
th = [th_; 2*pi-th_(end-1:-1:1)];

onesNth = ones(size(th_));
onesNr = ones(size(r_)); % not used much
eyeNth = spdiags(onesNth, 0, Nth, Nth);
eyeNr  = spdiags(onesNr', 0, Nr, Nr);
Nall = Nr*Nth;


% coef of H_rr: r^3
h_r2 = r_.^3;
h_rhat2 = onesNth * [r_.^3./r_rhat1.^2];

% coef of H_r: Pe(r-r^3)
h_r1 = Pe*(r_ - r_.^3);
h_rhat1 = onesNth * [h_r1./r_rhat1 - h_r2.*r_rhat2./r_rhat1.^3];

% coef of H_thth: r 
h_th2 = r_;
h_thhat2 = [1./th_thhat1.^2] * h_th2;

% coef of H_th: 2 Pe r sin(th)
%h_th1 = 2*Pe*sin(th_)*r_;
h_thhat1 = 2*Pe*[sin(th_)./th_thhat1]*r_ - [th_thhat2./th_thhat1.^3]*h_th2;

% coef of H: Pe r cos(th) + r/4 - Pe
h = (Pe*cos(th_) + 0.25)*r_ - Pe;

% old code
%h_r2 = R_.^3;
%h_r1 = Pe*(R_ - R_.^3);
%h_th2 = R_;
%h_th1 = 2*Pe*R_.*sin(TH_);
%h_0 = Pe*R_.*cos(TH_) + R_/4 - Pe;

h_rhat1(1,:) = 0; 
h_rhat2(1,:) = 0;
h_thhat1(1,:) = 1; 
h_thhat2(1,:) = 0;
h(1,:) = 0;

A = spdiags(h_rhat2(:),0,Nall,Nall) * kron(Drhat2, eyeNth) ...
    + spdiags(h_rhat1(:),0,Nall,Nall) * kron(Drhat1, eyeNth) ...
    + spdiags(h_thhat2(:),0,Nall,Nall) * kron(eyeNr, Dthhat2) ...
    + spdiags(h_thhat1(:),0,Nall,Nall) * kron(eyeNr, Dthhat1) ...
    + spdiags(h(:),0,Nall,Nall);

B = h_rhat2(:) .* kron(Drhat2_endcol, onesNth) ...
    + h_rhat1(:) .* kron(Drhat1_endcol, onesNth);

H_ = A\(-B);
H_ = reshape(H_,Nth,Nr);
H = [zeros(size(th_)) H_ onesNth];

[tmp,r_rhat1,tmp] = AdpmeshR(Pe, 1);
Sig = H*Drhat1_endrow'/r_rhat1 - 0.5; % Sig = H_rhat(r=1,th) - 1/2
Sig = Sig([1:end end-1:-1:1]); 

H = H([1:end end-1:-1:1],:); % th: [0 pi] --> [0 2pi]

C = H_ .* (onesNth*sqrt(1./r_)) .* exp( Pe*sin(th_/2).^2*(2-1./r_-r_) );
C = [zeros(size(th_)) C onesNth];
C = C([1:end end-1:-1:1],:); % th: [0 pi] --> [0 2pi]
