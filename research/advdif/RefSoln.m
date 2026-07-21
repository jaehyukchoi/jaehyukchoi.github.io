%
% RefSoln.m generates RefSoln.mat
%
%     RefSoln.mat stores the solutions for the references values of Pe 

%%%%
%%%%
%%%%

Nr = 72; Nth = 72;
% Make sure Nth is odd number to get \sigma at \th = \pi/2)

PeI = 10.^[-4:0.2:3];
n = length(PeI);

SigI = zeros(2*Nth-1,n);
HI = zeros(2*Nth-1,Nr+2,n);

for k = 1:n
  [SigI(:,k),r,th,T,HI(:,:,k)] = ConcSpectral(PeI(k),Nr,Nth);
  sprintf(' %7.3g%% done ...',k/n*100)
end

%%%%
%%%% Integration
%%%%
[Dth, th] = mycheb(2*Nth-2,[0 2*pi]);
Int = inv(Dth(2:end,2:end));

cdf_ = [zeros(1,n); Int*SigI(2:end,:)];
NuI = cdf_(end,:);
cdf_ = repop(cdf_,NuI,'/');

cdfI = th/(2*pi);
thI = zeros(size(cdf_));
for k = 1:n
  thI(:,k) = interp1(cdf_(:,k),th,cdfI,'spline');
end

save RefFlux r th Nr Nth PeI SigI NuI HI thI cdfI
