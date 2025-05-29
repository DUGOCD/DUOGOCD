function [f,varargout] = object_fun( x,f_num,x_num,fun )
% --------------------ZDT1--------------------
% [pop,~]=size(x);
if strcmp(fun,'ZDT1')
  f=[];
  f(1)=x(1);
  sum1=0;
  for i=2:x_num
      sum1 = sum1+x(i);
  end
  g=1+9*(sum1/(x_num-1));
  f(2)=g*(1-(f(1)/g)^0.5);
  %计算基准数值,真实PF值
            ff(:,1)    = (0:1/(100-1):1)';
            ff(:,2)    = 1-ff(:,1).^0.5;
            varargout = {ff};
end
% --------------------ZDT2--------------------
if strcmp(fun,'ZDT2')
  f=[];
  f(1)=x(1);
  sum1=0;
  for i=2:x_num
      sum1 = sum1+x(i);
  end
  g=1+9*(sum1/(x_num-1));
  f(2)=g*(1-(f(1)/g)^2);
   %计算基准数值,真实PF值
   ff(:,1)    = (0:1/(100-1):1)';
            ff(:,2)    = 1-ff(:,1).^2;
            varargout = {ff};
end
% --------------------ZDT3--------------------
if strcmp(fun,'ZDT3')
  f=[];
  f(1)=x(1);
  sum1=0;
  for i=2:x_num
      sum1 = sum1+x(i);
  end
  g=1+9*(sum1/(x_num-1));
  f(2)=g*(1-(f(1)/g)^0.5-(f(1)/g)*sin(10*pi*f(1)));
  %计算基准数值,真实PF值
   ff(:,1)    = (0:1/(100-1):1)';
            ff(:,2)    = 1-ff(:,1).^0.5-ff(:,1).*sin(10*pi*ff(:,1));
            varargout = {ff(NDSort(ff,1)==1,:)};
end
% --------------------ZDT4--------------------
if strcmp(fun,'ZDT4')
  f=[];
  f(1)=x(1);
  sum1=0;
  for i=2:x_num
      sum1 = sum1+(x(i)^2-10*cos(4*pi*x(i)));
  end
  g=1+9*10+sum1;
  f(2)=g*(1-(f(1)/g)^0.5);
   ff(:,1)    = (0:1/(100-1):1)';
            ff(:,2)    = 1-ff(:,1).^0.5;
            varargout = {ff};
end
% --------------------ZDT6--------------------
if strcmp(fun,'ZDT6')
  f=[];
  f(1)=1-(exp(-4*x(1)))*((sin(6*pi*x(1)))^6);
  sum1=0;
  for i=2:x_num
      sum1 = sum1+x(i);
  end
  g=1+9*((sum1/(x_num-1))^0.25);
  f(2)=g*(1-(f(1)/g)^2);
   minf1     = 0.280775;
            ff(:,1)    = (minf1:(1-minf1)/(100-1):1)';
            ff(:,2)    = 1-ff(:,1).^2;
            varargout = {ff};
end
% --------------------------------------------
% --------------------DTLZ1--------------------
if strcmp(fun,'DTLZ1')
  f=[];
   sum1=0;
  for i=3:x_num
      sum1 = sum1+((x(:,i)-0.5)^2-cos(20*pi*(x(:,i)-0.5)));
  end
  g=100*(x_num-2)+100*sum1;
%    g = 100*(5+sum1((x(:,f_num:end)-0.5).^2-cos(20.*pi.*(x(:,f_num:end)-0.5)),2));
  for i = 1 : f_num
                       f(:,i) = 0.5.*prod(x(:,1:f_num-i),2).*(1+g);
                        if i > 1
                            f(:,i) = f(:,i).*(1-x(:,f_num-i+1));
                        end
  end
  ff = UniformPoint(100,f_num)/2;
            varargout = {ff};
%  g      = 100*(x_num-f_num+1+sum1((x(:,f_num:end)-0.5).^2-cos(20.*pi.*(x(:,f_num:end)-0.5)),2));
%             f = 0.5*repmat(1+g,1,f_num).*fliplr(cumprod([ones(pop,1),x(:,1:f_num-1)],2)).*[ones(pop,1),1-x(:,f_num-1:-1:1)];
%   f(1)=(1+g)*x(1)*x(2);
%   f(2)=(1+g)*x(1)*(1-x(2));
%   f(3)=(1+g)*(1-x(1));
end
% --------------------------------------------
% --------------------DTLZ2--------------------
if strcmp(fun,'DTLZ2')
  f=[];
%   sum11=0;
%   for i=3:x_num
%       sum11 = sum11+(x(i))^2;
%   end
%   g=sum11;
%   f(1)=(1+g)*cos(x(1)*pi*0.5)*cos(x(2)*pi*0.5);
%   f(2)=(1+g)*cos(x(1)*pi*0.5)*sin(x(2)*pi*0.5);
%   f(3)=(1+g)*sin(x(1)*pi*0.5);
 g      = sum((x(:,f_num:end)-0.5).^2,2);
            f = repmat(1+g,1,f_num).*fliplr(cumprod([ones(size(g,1),1),cos(x(:,1:f_num-1)*pi/2)],2)).*[ones(size(g,1),1),sin(x(:,f_num-1:-1:1)*pi/2)];
   ff = UniformPoint(100,f_num);
            ff = ff./repmat(sqrt(sum(ff.^2,2)),1,f_num);
            % Return the reference points
            varargout = {ff};
end
if strcmp(fun,'DTLZ3')
     f=[];
     [N,~]=size(x);
            g = 100*(x_num-f_num+1+sum((x(:,f_num:end)-0.5).^2-cos(20.*pi.*(x(:,f_num:end)-0.5)),2));
            f = repmat(1+g,1,f_num).*fliplr(cumprod([ones(N,1),cos(x(:,1:f_num-1)*pi/2)],2)).*[ones(N,1),sin(x(:,f_num-1:-1:1)*pi/2)];           
            ff = UniformPoint(100,f_num);
            ff = ff./repmat(sqrt(sum(ff.^2,2)),1,f_num);
            varargout = {ff};
end
if strcmp(fun,'DTLZ4')
     f=[];
       x(:,1:f_num-1) = x(:,1:f_num-1).^100;
            g      = sum((x(:,f_num:end)-0.5).^2,2);
            f = repmat(1+g,1,f_num).*fliplr(cumprod([ones(size(g,1),1),cos(x(:,1:f_num-1)*pi/2)],2)).*[ones(size(g,1),1),sin(x(:,f_num-1:-1:1)*pi/2)];

            ff = UniformPoint(100,f_num);
            ff = ff./repmat(sqrt(sum(ff.^2,2)),1,f_num);
            varargout = {ff};
end
if strcmp(fun,'DTLZ5')
     f=[];
      g      = sum((x(:,f_num:end)-0.5).^2,2);
            Temp   = repmat(g,1,f_num-2);
            x(:,2:f_num-1) = (1+2*Temp.*x(:,2:f_num-1))./(2+2*Temp);
            f = repmat(1+g,1,f_num).*fliplr(cumprod([ones(size(g,1),1),cos(x(:,1:f_num-1)*pi/2)],2)).*[ones(size(g,1),1),sin(x(:,f_num-1:-1:1)*pi/2)];
            ff = [0:1/(100-1):1;1:-1/(100-1):0]';
            ff = ff./repmat(sqrt(sum(ff.^2,2)),1,size(ff,2));
            ff = [ff(:,ones(1,f_num-2)),ff];
            ff = ff./sqrt(2).^repmat([f_num-2,f_num-2:-1:0],size(ff,1),1);
            varargout  = {ff};
end
if strcmp(fun,'DTLZ6')
     f=[];
          g      = sum(x(:,f_num:end).^0.1,2);
            Temp   = repmat(g,1,f_num-2);
            x(:,2:f_num-1) = (1+2*Temp.*x(:,2:f_num-1))./(2+2*Temp);
            f = repmat(1+g,1,f_num).*fliplr(cumprod([ones(size(g,1),1),cos(x(:,1:f_num-1)*pi/2)],2)).*[ones(size(g,1),1),sin(x(:,f_num-1:-1:1)*pi/2)];
            ff = [0:1/(100-1):1;1:-1/(100-1):0]';
            ff = ff./repmat(sqrt(sum(ff.^2,2)),1,size(ff,2));
            ff = [ff(:,ones(1,f_num-2)),ff];
            ff = ff./sqrt(2).^repmat([f_num-2,f_num-2:-1:0],size(ff,1),1);
            varargout  = {ff};
end
if strcmp(fun,'DTLZ7')
     f=[]; 
     g               = 1+9*mean(x(:,f_num:end),2);
            f(:,1:f_num-1) = x(:,1:f_num-1);
            f(:,f_num)     = (1+g).*(f_num-sum(f(:,1:f_num-1)./(1+repmat(g,1,f_num-1)).*(1+sin(3*pi.*f(:,1:f_num-1))),2));

            interval     = [0,0.251412,0.631627,0.859401];
            median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
            X            = ReplicatePoint(100,f_num-1);
            X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
            X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
            ff            = [X,2*(f_num-sum(X/2.*(1+sin(3*pi.*X)),2))];
            varargout    = {ff};
end

if strcmp(fun,'LSMOP1')
      f=[]; 
     nk=5;
     M = f_num;
     PopDec = x;
     [N,D]  = size(PopDec);
     if isempty(M); M = 3; end
     if isempty(D); D = 100*M; end
%      D = 100*M
    c = 3.8*0.1*(1-0.1);
    for i = 1 : M-1
        c = [c,3.8.*c(end).*(1-c(end))];
    end
    sublen = floor(c./sum(c).*(D-M+1)/nk);
    len    = [0,cumsum(sublen*nk)];
%     [N,D]  = size(PopDec);
            PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
            G = zeros(N,M);
            for i = 1 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Sphere(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            for i = 2 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Sphere(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            G      = G./repmat(sublen,N,1)./nk;
            f = (1+G).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)]; 
            ff = UniformPoint(100,M);
            varargout = {ff};
end

if strcmp(fun,'LSMOP2')
       f=[]; 
     nk=5;
     M = f_num;
     PopDec = x;
     [N,D]  = size(PopDec);
     if isempty(M); M = 3; end
     if isempty(D); D = 100*M; end
     
    c = 3.8*0.1*(1-0.1);
    for i = 1 : M-1
        c = [c,3.8.*c(end).*(1-c(end))];
    end
    sublen = floor(c./sum(c).*(D-M+1)/nk);
    len    = [0,cumsum(sublen*nk)];
            
           PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
            G = zeros(N,M);
            for i = 1 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Griewank(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            for i = 2 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Schwefel(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            G      = G./repmat(sublen,N,1)./nk;
            f = (1+G).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)]; 
            ff = UniformPoint(100,M);
            varargout = {ff};
end
if strcmp(fun,'LSMOP3')
       f=[]; 
     nk=5;
     M = f_num;
     PopDec = x;
     [N,D]  = size(PopDec);
     if isempty(M); M = 3; end
     if isempty(D); D = 100*M; end
     
    c = 3.8*0.1*(1-0.1);
    for i = 1 : M-1
        c = [c,3.8.*c(end).*(1-c(end))];
    end
    sublen = floor(c./sum(c).*(D-M+1)/nk);
    len    = [0,cumsum(sublen*nk)];
            
            
           PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
            G = zeros(N,M);
            for i = 1 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Rastrigin(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            for i = 2 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Rosenbrock(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            G      = G./repmat(sublen,N,1)./nk;
            f = (1+G).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)]; 
            ff = UniformPoint(100,M);
            varargout = {ff};
end
if strcmp(fun,'LSMOP4')
      f=[]; 
     nk=5;
     M = f_num;
     PopDec = x;
     [N,D]  = size(PopDec);
     if isempty(M); M = 3; end
     if isempty(D); D = 100*M; end
     
    c = 3.8*0.1*(1-0.1);
    for i = 1 : M-1
        c = [c,3.8.*c(end).*(1-c(end))];
    end
    sublen = floor(c./sum(c).*(D-M+1)/nk);
    len    = [0,cumsum(sublen*nk)];
           PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
            G = zeros(N,M);
           for i = 1 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Ackley(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            for i = 2 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Griewank(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            G      = G./repmat(sublen,N,1)./nk;
            f = (1+G).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)]; 
            ff = UniformPoint(100,M);
            varargout = {ff};
end
if strcmp(fun,'LSMOP5')
       f=[]; 
     nk=5;
     M = f_num;
     PopDec = x;
     [N,D]  = size(PopDec);
     if isempty(M); M = 3; end
     if isempty(D); D = 100*M; end
     
    c = 3.8*0.1*(1-0.1);
    for i = 1 : M-1
        c = [c,3.8.*c(end).*(1-c(end))];
    end
    sublen = floor(c./sum(c).*(D-M+1)/nk);
    len    = [0,cumsum(sublen*nk)];
           PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
            G = zeros(N,M);
           for i = 1 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Sphere(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            for i = 2 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Sphere(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            G      = G./repmat(sublen,N,1)./nk;
            f = (1+G+[G(:,2:end),zeros(N,1)]).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
            ff = UniformPoint(100,M);
            ff=ff./repmat(sqrt(sum(ff.^2,2)),1,M);
            varargout = {ff};
end
if strcmp(fun,'LSMOP6')
       f=[]; 
     nk=5;
     M = f_num;
     PopDec = x;
     [N,D]  = size(PopDec);
     if isempty(M); M = 3; end
     if isempty(D); D = 100*M; end
     
    c = 3.8*0.1*(1-0.1);
    for i = 1 : M-1
        c = [c,3.8.*c(end).*(1-c(end))];
    end
    sublen = floor(c./sum(c).*(D-M+1)/nk);
    len    = [0,cumsum(sublen*nk)];
            
           PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
            G = zeros(N,M);
           for i = 1 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Rosenbrock(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            for i = 2 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Schwefel(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            G      = G./repmat(sublen,N,1)./nk;
            f = (1+G+[G(:,2:end),zeros(N,1)]).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
            ff = UniformPoint(100,M);
            ff=ff./repmat(sqrt(sum(ff.^2,2)),1,M);
            varargout = {ff};
end
if strcmp(fun,'LSMOP7')
       f=[]; 
     nk=5;
     M = f_num;
     PopDec = x;
     [N,D]  = size(PopDec);
     if isempty(M); M = 3; end
     if isempty(D); D = 100*M; end
     
    c = 3.8*0.1*(1-0.1);
    for i = 1 : M-1
        c = [c,3.8.*c(end).*(1-c(end))];
    end
    sublen = floor(c./sum(c).*(D-M+1)/nk);
    len    = [0,cumsum(sublen*nk)];
            
           PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
            G = zeros(N,M);
           for i = 1 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Ackley(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            for i = 2 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Rosenbrock(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            G      = G./repmat(sublen,N,1)./nk;
            f = (1+G+[G(:,2:end),zeros(N,1)]).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
            
            ff = UniformPoint(100,M);
            ff=ff./repmat(sqrt(sum(ff.^2,2)),1,M);
            varargout = {ff};
end
if strcmp(fun,'LSMOP8')
     f=[]; 
     nk=5;
     M = f_num;
     PopDec = x;
     [N,D]  = size(PopDec);
     if isempty(M); M = 3; end
     if isempty(D); D = 100*M; end
     
    c = 3.8*0.1*(1-0.1);
    for i = 1 : M-1
        c = [c,3.8.*c(end).*(1-c(end))];
    end
    sublen = floor(c./sum(c).*(D-M+1)/nk);
    len    = [0,cumsum(sublen*nk)];
            
           PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
            G = zeros(N,M);
           for i = 1 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Griewank(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            for i = 2 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Sphere(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            G      = G./repmat(sublen,N,1)./nk;
            f = (1+G+[G(:,2:end),zeros(N,1)]).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
            ff = UniformPoint(100,M);
            ff=ff./repmat(sqrt(sum(ff.^2,2)),1,M);
            varargout = {ff};
end
if strcmp(fun,'LSMOP9')
      f=[]; 
     nk=5;
     M = f_num;
     PopDec = x;
     [N,D]  = size(PopDec);
     if isempty(M); M = 3; end
     if isempty(D); D = 100*M; end
%      D = 100*M
    c = 3.8*0.1*(1-0.1);
    for i = 1 : M-1
        c = [c,3.8.*c(end).*(1-c(end))];
    end
%      sublen = ceil(round(c./sum(c).*D)/nk);
      sublen = floor(c./sum(c).*(D-M+1)/nk);
    len    = [0,cumsum(sublen*nk)];
%      PopDec = x;
%             [N,D]  = size(PopDec);
            
           PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
            G = zeros(N,M);
            for i = 1 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Sphere(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            for i = 2 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Ackley(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            G = 1 + sum(G./repmat(sublen,N,1)./nk,2);
            f(:,1:M-1) = PopDec(:,1:M-1);
            f(:,M)     = (1+G).*(M-sum(f(:,1:M-1)./(1+repmat(G,1,M-1)).*(1+sin(3*pi.*f(:,1:M-1))),2));
             interval     = [0,0.251412,0.631627,0.859401];
            median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
            X            = ReplicatePoint(100,f_num-1);
            X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
            X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
            ff            = [X,2*(f_num-sum(X/2.*(1+sin(3*pi.*X)),2))];
            varargout    = {ff};
end


if strcmp(fun,'WFG1')
      f=[]; 
     PopDec = x;
            [N,D]  = size(PopDec);
            M      = f_num;
            K = 2*(M-1);
            L = D - K;
            D = 1;
            S = 2 : 2 : 2*M;
            A = ones(1,M-1);

            z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);

            t1 = zeros(N,K+L);
            t1(:,1:K)     = z01(:,1:K);
            t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);

            t2 = zeros(N,K+L);
            t2(:,1:K)     = t1(:,1:K);
            t2(:,K+1:end) = b_flat(t1(:,K+1:end),0.8,0.75,0.85);

            t3 = zeros(N,K+L);
            t3 = b_poly(t2,0.02);

            t4 = zeros(N,M);
            for i = 1 : M-1
                t4(:,i) = r_sum(t3(:,(i-1)*K/(M-1)+1:i*K/(M-1)),2*((i-1)*K/(M-1)+1):2:2*i*K/(M-1));
            end
            t4(:,M) = r_sum(t3(:,K+1:K+L),2*(K+1):2:2*(K+L));

            x = zeros(N,M);
            for i = 1 : M-1
                x(:,i) = max(t4(:,M),A(i)).*(t4(:,i)-0.5)+0.5;
            end
            x(:,M) = t4(:,M);

            h      = convex(x);
            h(:,M) = mixed(x);
            f = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
            
            h = UniformPoint(100,M);
            for i = 1 : size(h,1)
                c = ones(1,M);
                k = find(h(i,:)~=0,1);
                for j = k+1 : M
                    temp     = h(i,j)/h(i,k)*prod(1-c(M-j+2:M-k));
                    c(M-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
                end
                for j = 1 : M
                    h(i,j) = prod(1-c(1:M-j)).*(1-sqrt(1-c(M-j+1)^2));
                end
                temp   = acos(c(1))*2/pi;                   
                h(i,M) = 1 - temp - cos(10*pi*temp+pi/2)/10/pi;
            end
            h = repmat(2:2:2*M,size(h,1),1).*h;
            varargout = {h};
end
% --------------------------------------------
end
function W = ReplicatePoint(SampleNum,M)
    if M > 1
        SampleNum = (ceil(SampleNum^(1/M)))^M;
        Gap       = 0:1/(SampleNum^(1/M)-1):1;
        eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
        eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
    else
        W = (0:1/(SampleNum-1):1)';
    end
end
function f = Sphere(x)
    f = sum(x.^2,2);
end
function f = Griewank(x)
    f = sum(x.^2,2)./4000 - prod(cos(x./repmat(sqrt(1:size(x,2)),size(x,1),1)),2) + 1;
end

function f = Schwefel(x)
    f = max(abs(x),[],2);
end
function f = Rastrigin(x)
    f = sum(x.^2-10.*cos(2.*pi.*x)+10,2);
end

function f = Rosenbrock(x)
    f = sum(100.*(x(:,1:size(x,2)-1).^2-x(:,2:size(x,2))).^2+(x(:,1:size(x,2)-1)-1).^2,2);
end
function f = Ackley(x)
    f = 20-20.*exp(-0.2.*sqrt(sum(x.^2,2)./size(x,2)))-exp(sum(cos(2.*pi.*x),2)./size(x,2))+exp(1);
end


function Output = s_linear(y,A)
    Output = abs(y-A)./abs(floor(A-y)+A);
end

function Output = b_flat(y,A,B,C)
    Output = A+min(0,floor(y-B))*A.*(B-y)/B-min(0,floor(C-y))*(1-A).*(y-C)/(1-C);
    Output = roundn(Output,-6);
end

function Output = b_poly(y,a)
    Output = y.^a;
end

function Output = r_sum(y,w)
    Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = convex(x)
    Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
end

function Output = mixed(x)
    Output = 1-x(:,1)-cos(10*pi*x(:,1)+pi/2)/10/pi;
end
