Pe = 1;

Clevel = linspace(0,1+1e-5,15);

[C, r, th] = GetRefSoln(Pe);

% frame
rr = cos([0:72]*pi/72)/2 + 1/2;
figure, plot([rr-i*0.2*pi; rr+i*2.2*pi],'k','linewidth',1)
hold on
plot([i*2*pi*rr-0.1; + i*2*pi*rr+1.1],'k','linewidth',1)

contourf(r,th,C,Clevel), axis square
axis([-0.1 1.1 -0.2*pi 2.2*pi]);
set(gca,'ytick',[0 pi 2*pi],'yticklabel',[{'0'}, {'1pi'}, {'2pi'}])
xlabel('r','fontsize',24);
ylabel('\theta','fontsize',24);

