Pe = 1;

Clevel = linspace(0,1+1e-5,15);
Flevel = setdiff([-4:1:4],0);

[C_, r_, th] = GetRefSoln(Pe);
r = r_(2:end);
C = C_(:,2:end);

Zin = exp(i*th) * r;
Zpl = Zin + 1./Zin;

figure, contourf(real(Zin),imag(Zin),C,Clevel), axis equal
hold on
[c,h]=contour(real(Zin),imag(Zin),imag(Zpl), Flevel,'w');
set(h,'linewidth',2);
plot(complex([-1 1]'),'w','linewidth',2);
axis([-1.1 1.1 -1.1 1.1]);
xlabel('Re z','fontsize',24)
ylabel('Im z','fontsize',24)
