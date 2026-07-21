function [C, varargout] = GetRefSoln(Pe,rI,thI);
%GETREFSOLN gives the solution from the interpolation of REFSOLN.mat
%
%     C = GetRefSoln(Pe,r,th);
%     [C, Z] = GetRefSoln(Pe);
%     [C, r, th] = GetRefSoln(Pe);

if nargin < 1
  error('We need Pe as the first argument at least.')
end

load RefSoln r th HI

k = nnz(log(PeI) < log(Pe));
if k == length(PeI)
  H = HI(:,:,end);
else
  a1 = (LogPeI(k+1)-LogPe)/(LogPeI(k+1)-LogPeI(k));
  a2 = (LogPe-LogPeI(k))/(LogPeI(k+1)-LogPeI(k));
  H = a1*HI(:,:,k) + a2*HI(:,:,k+1);
end

if ~exist('rI') | isempty(rI), rI = r; end
if ~exist('thI') | isempty(thI), thI = th; end

rI_ = rI;
I = (rI>1); rI_(I) = 1./rI(I);

if nargin > 1
  H = interp2(r, th, H, rI_, thI, 'spline');
end

if nargout == 2
  varargout{1} = exp(i*thI)*rI;
elseif nargout == 3
  varargout{1} = rI; varargout{2} = thI; 
end

I = (rI_>0);
C = zeros(length(thI),length(rI_));
C(:,I) = H(:,I) .* [ones(size(thI))*sqrt(1./rI_(I))] ... 
	 .* exp( Pe*sin(thI/2).^2 * (2-rI_(I)-1./rI_(I)));
