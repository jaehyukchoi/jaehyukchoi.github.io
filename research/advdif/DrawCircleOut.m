Pe = 1;

AxisLimit = [-3 8 -4.5 4.5];
Clevel = linspace(0,1+1e-5,15);
Flevel = setdiff([-4:1:4],0);

[C_, r_, th] = GetRefSoln(Pe);
r = r_(2:end);
C = C_(:,2:end);

Zout = exp(i*th) * (1./r);
Zpl = Zout + 1./Zout;

figure, contourf(real(Zout), imag(Zout), C, Clevel), axis equal
axis(AxisLimit);
set(gca,'xtick',[-4:2:6],'ytick',[-4:2:4]);
hold on
[c,h] = contour(real(Zout), imag(Zout), imag(Zpl), Flevel,'w');
set(h,'linewidth',2);
plot(complex([-10 -1; 1 10]'),'w','linewidth',2);
