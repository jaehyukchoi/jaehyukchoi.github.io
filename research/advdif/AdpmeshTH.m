function [th,th_psi1,th_psi2] = AdpmeshTH(Pe,psi)
% ADPMESHTH generates adaptive mesh for radial direction.
% 
%     [th,th_psi1,th_psi2] = adpmesh_th(Pe, th)
%
%     th(psi)   = x - (1-Pe) sin(psi)
%     th'(psi)  = 1 - (1-Pe) cos(psi)
%     th''(psi) = (1-Pe) sin(psi)
%
%     ADPMESHR and ADPMESHTH are used in FLUXADPMESH

HIGH = 10;

if HIGH < Pe
  a = 1-1/sqrt(Pe);
  th = psi - a*sin(psi);
  th_psi1 = 1 - a*cos(psi);
  th_psi2 = a*sin(psi);
else
  th = psi;
  th_psi1 = ones(size(psi));
  th_psi2 = zeros(size(psi));
end
