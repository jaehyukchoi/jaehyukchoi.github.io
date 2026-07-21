function Nu = Nussel(Pe)

PE_LOW = 0.001;

Nu = zeros(size(Pe));

[th, I] = myclencurt(200, [0 2*pi]);

for k = 1:length(Pe)
  if Pe(k) > PE_LOW
    Sig = FluxAsymH(Pe(k),th');
  else
    Sig = FluxAsymL(Pe(k),th');
  end
  Nu(k) = Sig' * I;
end
