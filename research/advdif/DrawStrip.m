Pe = 1;

AxisLimit = [-3 8 -4.5 4.5];
Clevel = linspace(0,1+1e-5,15);
Flevel = setdiff([-4:1:4],0);

[C_, r_, th] = GetRefSoln(Pe);
r = r_(2:end);
C = C_(:,2:end);

Zout = exp(i*th) * (1./r);
Zpl = Zout + 1./Zout;

figure, contourf(real(Zpl)/2,imag(Zpl)/2,C,Clevel), axis equal
hold on
[c,h]=contour(real(Zpl)/2,imag(Zpl)/2,imag(Zpl),Flevel,'w');
set(h,'linewidth',2);

plot(complex([-10 -1; 1 20]'),'w','linewidth',2)
plot(complex([-1 1]'),'linewidth',4,'color',[1 1 1])
set(gca,'xtick',[-1 0 1],'ytick',[])
axis(AxisLimit/2);

