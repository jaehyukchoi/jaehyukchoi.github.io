Pe = [1];
AxisLimit = [-3 8 -4.5 4.5]; % JFM 2004
%AxisLimit = [-2.5 4 -2.25 2.25]; % PRL 2003
Gray = [1 1 1]*0.35;
Clevel = linspace(0,1+1e-5,15); %Make sure to inlude 1
Flevel = setdiff([-4:1:4],0); %setdiff([-4:1:4],0)

%%%% Which to plot
Out = 0;
In = 0;
InOut = 0;
Tail = 0;
Plate = 0;
Polar = 1;
Polygon = 0;

for m = 1:length(Pe)

  [C_, r_, th] = GetRefSoln(Pe(m));
  r = r_(2:end);
  C = C_(:,2:end);
  
  Zin = exp(i*th) * r;
  Zin_ = exp(i*th) * r_;
  Zout = exp(i*th) * (1./r);
  Zpl = Zin + 1./Zin;
  
  %%%% Inside Circle
  if In
    figure, zcontourf(Zin_,C_,Clevel), axis equal
    hold on
    [c,h]=zcontour(Zin,Zpl,Flevel,'k');
    set(h,'linewidth',2);
    plot(complex([-1 1]'),'k','linewidth',2);
    axis([-1.1 1.1 -1.1 1.1]);
    xlabel('Re z','fontsize',24)
    ylabel('Im z','fontsize',24)
  end

  %%%% Outside
  if Out
    figure, zcontourf(Zout,C,Clevel), axis equal
    axis(AxisLimit);
    set(gca,'xtick',[-4:2:6],'ytick',[-4:2:4]);
    hold on
    [c,h]=zcontour(Zout,Zpl, Flevel,'k');
    set(h,'linewidth',2);
    plot(complex([-10 -1; 1 10]'),'k','linewidth',2);
  end

  %%%% Outside Circle: Tail
  if Tail
    figure, zcontourf(Zout,C,Clevel), axis equal
    hold on
    [c,h]=zcontour(Zout,Zpl, Flevel,'k');
    set(h,'linewidth',2);
    plot(complex([-10 -1; 1 10]'),'k','linewidth',2);
    set(gca,'xtick',[],'ytick',[]);
    axis([-3 10 -1.8 1.8]);
  end
  
  %%%% Finite Plate
  if Plate
    figure, zcontourf(Zpl/2,C,Clevel), axis equal
    hold on
    [c,h]=zcontour(Zpl/2,Zpl,Flevel,'k');
    set(h,'linewidth',2);
    
    plot(complex([-10 -1; 1 20]'),'k','linewidth',2)
    plot(complex([-1 1]'),'linewidth',4,'color',[1 1 1])
    set(gca,'xtick',[-1 0 1],'ytick',[])
    axis(AxisLimit/2);
    %colorbar('horiz')
  end
  
  %%%% Polygon using SC maps
  if Polygon
    load SCmaps f1
    SCmap = f1;
    Zpol = feval(SCmap,Zin);
    
    figure, zcontourf(Zpol,C,Clevel), axis equal
    hold on
    [c,h]=zcontour(Zpol,Zpl,Flevel,'k');
    set(h,'linewidth',2);
    
    plot(complex([-10 feval(SCmap,-1); feval(SCmap,1) 20]'),'k','linewidth',2)
    plot(complex([-1 1]'),'linewidth',4,'color',[1 1 1])
    set(gca,'xtick',[-1 0 1],'ytick',[])
    axis(AxisLimit*0.83462684);
    %colorbar('horiz')
  end

  %%%% Polar Coordinate
  if Polar
    rr = cos([0:pi/20:pi])/2 + 1/2;
    
    figure, plot([rr-i*0.2*pi; rr+i*2.2*pi],'k','linewidth',1)
    hold on
    plot([i*2*pi*rr-0.1; + i*2*pi*rr+1.1],'k','linewidth',1)
    contourf(r_,th,C_,Clevel), axis square
    %[c,h]= contour(r,th,real(Zpl), [0 Flevel],'k');
    %set(h,'linewidth',2);
    axis([-0.1 1.1 -0.2*pi 2.2*pi]);
    set(gca,'ytick',[0 pi 2*pi],'yticklabel',[{'0'}, {'1pi'}, {'2pi'}])
    xlabel('r','fontsize',24);
    ylabel('\theta','fontsize',24);
  end
end
