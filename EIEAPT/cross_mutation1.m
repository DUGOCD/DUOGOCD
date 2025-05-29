function chromo_offspring = cross_mutation1( chromo_parent,x_num,x_min,x_max,pc,pm,yita1,yita2 )
    ParentDec = chromo_parent(:,1:x_num);
    [N,D]     = size(ParentDec);

    %% Simulated binary crossover
    
    
    Parent1Dec = ParentDec(1:N/2,:);
    Parent2Dec = ParentDec(N/2+1:end,:);
    beta = zeros(N/2,D);
    mu   = rand(N/2,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(yita1+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(yita1+1));
    beta = beta.*(-1).^randi([0,1],N/2,D);
    beta(rand(N/2,D)<0.5) = 1;
    beta(repmat(rand(N/2,1)>pc,1,D)) = 1;
    OffspringDec = [(Parent1Dec+Parent2Dec)/2+beta.*(Parent1Dec-Parent2Dec)/2
                    (Parent1Dec+Parent2Dec)/2-beta.*(Parent1Dec-Parent2Dec)/2];

    %% Polynomial mutation
    Lower = repmat(x_min,N,1);
    Upper = repmat(x_max,N,1);
    Site  = rand(N,D) <pm;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                         (1-(OffspringDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(yita2+1)).^(1/(yita2+1))-1);
    temp = Site & mu>0.5; 
    OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                         (1-(Upper(temp)-OffspringDec(temp))./(Upper(temp)-Lower(temp))).^(yita2+1)).^(1/(yita2+1)));
       
%     OffspringDec( OffspringDec(:,1)<0)=0;
%     OffspringDec( OffspringDec(:,1)>1)=1;
%     for i=2:x_num
%     OffspringDec( OffspringDec(:,i)<-5)=-5;
%     OffspringDec( OffspringDec(:,i)>5)=5;
%     end
for i=1:N
  for j=1:x_num
    if(OffspringDec(i,j)>x_max(j))
               OffspringDec(i,j)=x_max(j);
           elseif(OffspringDec(i,j)<x_min(j))
               OffspringDec(i,j)=x_min(j);
    end
  end
end
%     for i=1:N
%         [OffspringDec(i,(x_num+1):(x_num+f_num)),~]=object_fun(OffspringDec(i,1:x_num),f_num,x_num,fun);
%     end
% %ģ������ƽ��������ʽ����
% [pop,~]=size(chromo_parent);
% suoyin=1;
% for i=1:pop
%    %%%ģ������ƽ���
%    %��ʼ���Ӵ���Ⱥ
%    %���ѡȡ������������
%    SUJI=randperm(length(chromo_parent(:,1)),2)
%    parent_1=SUJI(1);
%    parent_2=SUJI(2);
% %    parent_1=round(pop*rand(1));
% %    if (parent_1<1)
% %        parent_1=1;
% %    end
% %    parent_2=round(pop*rand(1));
% %    if (parent_2<1)
% %        parent_2=1;
% %    end
% %    %ȷ�������������岻��ͬһ��
% %    while isequal(chromo_parent(parent_1,:),chromo_parent(parent_2,:))
% %        parent_2=round(pop*rand(1));
% %        if(parent_2<1)
% %            parent_2=1;
% %        end
% %    end
%    chromo_parent_1=chromo_parent(parent_1,:);
%    chromo_parent_2=chromo_parent(parent_2,:);
%    off_1=chromo_parent_1;
%    off_2=chromo_parent_2;%���ǲ��Ǿ����������أ��ðɣ����Ѿ��ĳ���2
%    if(rand(1)<pc)
%        %����ģ������ƽ���
%        u1=zeros(1,x_num);
%        gama=zeros(1,x_num);
%        for j=1:x_num
%            u1(j)=rand(1);
%            if u1(j)<0.5
%                gama(j)=(2*u1(j))^(1/(yita1+1));
%            else
%                gama(j)=(1/(2*(1-u1(j))))^(1/(yita1+1));
%            end
%            off_1(j)=0.5*((1+gama(j))*chromo_parent_1(j)+(1-gama(j))*chromo_parent_2(j));
%            off_2(j)=0.5*((1-gama(j))*chromo_parent_1(j)+(1+gama(j))*chromo_parent_2(j));
%            %ʹ�Ӵ��ڶ�������
%            if(off_1(j)>x_max(j))
%                off_1(j)=x_max(j);
%            elseif(off_1(j)<x_min(j))
%                off_1(j)=x_min(j);
%            end
%            if(off_2(j)>x_max(j))
%                off_2(j)=x_max(j);
%            elseif(off_2(j)<x_min(j))
%                off_2(j)=x_min(j);
%            end
%        end
%        %�����Ӵ������Ŀ�꺯��ֵ
%        [off_1(1,(x_num+1):(x_num+f_num)),~]=object_fun(off_1,f_num,x_num,fun);
%        [off_2(1,(x_num+1):(x_num+f_num)),~]=object_fun(off_2,f_num,x_num,fun);
%    end
%    %%%����ʽ����
%    if(rand(1)<pm)
%        u2=zeros(1,x_num);
%        delta=zeros(1,x_num);
%        for j=1:x_num
%            u2(j)=rand(1);
%            if(u2(j)<0.5)
%                delta(j)=(2*u2(j))^(1/(yita2+1))-1;
%            else
%                delta(j)=1-(2*(1-u2(j)))^(1/(yita2+1));
%            end
%            off_1(j)=off_1(j)+delta(j);
%            %ʹ�Ӵ��ڶ�������
%            if(off_1(j)>x_max(j))
%                off_1(j)=x_max(j);
%            elseif(off_1(j)<x_min(j))
%                off_1(j)=x_min(j);
%            end
%        end
%        %�����Ӵ������Ŀ�꺯��ֵ
%        [off_1(1,(x_num+1):(x_num+f_num)),~]=object_fun(off_1,f_num,x_num,fun);
%    end
%    if(rand(1)<pm)
%        u2=zeros(1,x_num);
%        delta=zeros(1,x_num);
%        for j=1:x_num
%            u2(j)=rand(1);
%            if(u2(j)<0.5)
%                delta(j)=(2*u2(j))^(1/(yita2+1))-1;
%            else
%                delta(j)=1-(2*(1-u2(j)))^(1/(yita2+1));
%            end
%            off_2(j)=off_2(j)+delta(j);
%            %ʹ�Ӵ��ڶ�������
%            if(off_2(j)>x_max(j))
%                off_2(j)=x_max(j);
%            elseif(off_2(j)<x_min(j))
%                off_2(j)=x_min(j);
%            end
%        end
%        %�����Ӵ������Ŀ�꺯��ֵ
%        [off_2(1,(x_num+1):(x_num+f_num)),~]=object_fun(off_2,f_num,x_num,fun);
%    end
%    off(suoyin,:)=off_1;
%    off(suoyin+1,:)=off_2;
%    suoyin=suoyin+2;
% end
% pick = randperm(2*pop, pop)'; % �� 2*NP �����ѡ�� NP ��
%    chromo_offspring= off(pick, :);
OffspringDec(:,x_num+1)=chromo_parent(:,x_num+1);
 chromo_offspring=OffspringDec;
end
