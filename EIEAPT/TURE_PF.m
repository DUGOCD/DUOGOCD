function TURE_PF(fun,PF)
if strcmp(fun,'DTLZ1')
a = linspace(0,1,10)';
 R = {a*a'/2,a*(1-a')/2,(1-a)*ones(size(a'))/2};
  h=surf(R{1},R{2},R{3});
              set(h,'facecolor','none','edgecolor',[.7 .7 .7])
end
if strcmp(fun,'DTLZ2')
     a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
                 h=surf(R{1},R{2},R{3});
              set(h,'facecolor','none','edgecolor',[.7 .7 .7])
end
if strcmp(fun,'DTLZ3')
     a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
                 h=surf(R{1},R{2},R{3});
              set(h,'facecolor','none','edgecolor',[.7 .7 .7])
end
if strcmp(fun,'DTLZ4')
    a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
                 h=surf(R{1},R{2},R{3});
              set(h,'facecolor','none','edgecolor',[.7 .7 .7])
end
if strcmp(fun,'DTLZ5')
    plot3(PF(:,1),PF(:,2),PF(:,3),'k-');
end
if strcmp(fun,'DTLZ6')
   plot3(PF(:,1),PF(:,2),PF(:,3),'k-');
end
if strcmp(fun,'DTLZ7')
    [x,y]  = meshgrid(linspace(0,1,20));
                z      = 2*(3-x/2.*(1+sin(3*pi*x))-y/2.*(1+sin(3*pi*y)));
                nd     = reshape(NDSort([x(:),y(:),z(:)],1)==1,size(z));
                z(~nd) = nan;
                R      = {x,y,z}; 
                 h=surf(R{1},R{2},R{3});
              set(h,'facecolor','none','edgecolor',[.7 .7 .7])
end
if strcmp(fun,'LSMOP1')
     a = linspace(0,1,10)';
                R = {a*a',a*(1-a'),(1-a)*ones(size(a'))};
     h=surf(R{1},R{2},R{3});
              set(h,'facecolor','none','edgecolor',[.7 .7 .7])
end
if strcmp(fun,'LSMOP2')
    a = linspace(0,1,10)';
                R = {a*a',a*(1-a'),(1-a)*ones(size(a'))};
     h=surf(R{1},R{2},R{3});
              set(h,'facecolor','none','edgecolor',[.7 .7 .7])
end
if strcmp(fun,'LSMOP3')
    a = linspace(0,1,10)';
                R = {a*a',a*(1-a'),(1-a)*ones(size(a'))};
     h=surf(R{1},R{2},R{3});
              set(h,'facecolor','none','edgecolor',[.7 .7 .7])
end
if strcmp(fun,'LSMOP4')
    a = linspace(0,1,10)';
                R = {a*a',a*(1-a'),(1-a)*ones(size(a'))};
     h=surf(R{1},R{2},R{3});
              set(h,'facecolor','none','edgecolor',[.7 .7 .7])
end
if strcmp(fun,'LSMOP5')
      a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
     h=surf(R{1},R{2},R{3});
              set(h,'facecolor','none','edgecolor',[.7 .7 .7])
end
if strcmp(fun,'LSMOP6')
      a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
     h=surf(R{1},R{2},R{3});
              set(h,'facecolor','none','edgecolor',[.7 .7 .7])
end
if strcmp(fun,'LSMOP7')
      a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
     h=surf(R{1},R{2},R{3});
              set(h,'facecolor','none','edgecolor',[.7 .7 .7])
end
if strcmp(fun,'LSMOP8')
     a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
     h=surf(R{1},R{2},R{3});
              set(h,'facecolor','none','edgecolor',[.7 .7 .7])
end
if strcmp(fun,'LSMOP9')
     [x,y]  = meshgrid(linspace(0,1,20));
                z      = 2*(3-x/2.*(1+sin(3*pi*x))-y/2.*(1+sin(3*pi*y)));
                nd     = reshape(NDSort([x(:),y(:),z(:)],1)==1,size(z));
                z(~nd) = nan;
                R      = {x,y,z};
     h=surf(R{1},R{2},R{3});
              set(h,'facecolor','none','edgecolor',[.7 .7 .7])
end

end