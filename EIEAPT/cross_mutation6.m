function chromo_offspring = cross_mutation6( chromo_parent,f_num,x_num,x_min,x_max,pc,pm,yita1,yita2,fun,F1 )
   %% 坏的向好的学，好的互相学，加上了双向向量,并且进行了改进,四象之力
    ParentDec = chromo_parent(:,1:x_num);
    [N,D]     = size(ParentDec);
     numbers = 1:N;
     complementNumbers = setdiff(numbers, F1(1).ss);
     if ~isempty(complementNumbers)
     FN=ParentDec(F1(1).ss,:);
     FD=ParentDec(complementNumbers,:);
      if size(F1(1).ss,2)>=2
     for i=1:size(F1(1).ss,2)
         
         distances = pdist2(FN(i,:), FN);
         SOR=sort(distances);
         SORS(i)=SOR(2);
     end
     DAV=sum(SORS)/size(F1(1).ss,2);
    else
        SORS=1
     DAV=2;
    end
     for i=1:size(F1(1).ss,2)
         if SORS(i)>DAV
           if size(F1(1).ss,2)==1
             if rand<0.5
             DD=(FN(i,:)-x_max)/norm(FN(i,:)-x_max);
             else
             DD=(FN(i,:)-x_min)/norm(FN(i,:)-x_min);
             end            
          gama=  randn(1, 1);
          OffspringDec(i,:)=FN(i,:)+DD*gama;
           else
            NUM=1:size(F1(1).ss,2);
           HBNUM = setdiff(NUM, i);
           numElements = length(HBNUM);
          % 生成一个随机索引，范围从 1 到数字集合的大小
          randomIndex = randi(numElements);
          selectedNumber = HBNUM(randomIndex);
          DD=(FN(i,:)-FN(selectedNumber,:))/norm(FN(i,:)-FN(selectedNumber,:));
          for k=1:size(F1(1).ss,2)
             namda(k)=(FN(i,:)-FN(k,:))*DD';      
          end
         variance = var(namda);
         std_dev = sqrt(variance);
         gama= std_dev * randn(1, 1);
         OffspringDec(i,:)=FN(i,:)+DD*gama;
           end
         else
              if rand<0.5
             DD=(FN(i,:)-x_max)/norm(FN(i,:)-x_max);
             else
             DD=(FN(i,:)-x_min)/norm(FN(i,:)-x_min);
             end
           if  size(F1(1).ss,2)==1
               for k=1:size(complementNumbers,2)
                namda(k)=(FN(1,:)-FD(k,:))*DD';      
               end
           else
          for k=1:size(F1(1).ss,2)
             namda(k)=(FN(i,:)-FN(k,:))*DD';      
          end
           end
         variance = var(namda);
         std_dev = sqrt(variance);
         gama= std_dev * randn(1, 1);
%          gama=abs(gama);
         OffspringDec(i,:)=FN(i,:)+DD*gama;
         end
     end 
      for i=1:size(complementNumbers,2)
    %%  找到方向
        % 生成一个随机索引，范围从 1 到数字集合的大小
         selectedNumber = randi([1, size(F1(1).ss,2)]);
%          DD=(FN(selectedNumber,:)-FD(i,:))/norm(FN(selectedNumber,:)-FD(i,:));
          kkkk=rand;
         if kkkk<0.5
             v=FD(i,:)-x_max;
             null_space=null(v);
             CCC=randi([1, size(null_space,2)]);
             DD=null_space(:,1)'/norm(null_space(:,1)');%              DD=(FN(i,:)-x_max)/norm(FN(i,:)-x_max);
         elseif    0.5<= kkkk<1
              v=FD(i,:)-x_min;
             null_space=null(v);
%              ks=null(kkk,'r');
             CCC=randi([1, size(null_space,2)]);
             DD=null_space(:,1)'/norm(null_space(:,1)');%              DD=(FN(i,:)-x_max)/norm(FN(i,:)-x_max);
         else
             DD=(FN(selectedNumber,:)-FD(i,:))/norm(FN(selectedNumber,:)-FD(i,:));
         end
         for k=1:size(complementNumbers,2)
            namda(k)=(FN(selectedNumber,:)-FD(k,:))*DD';      
         end
         variance = var(namda);
         std_dev = sqrt(variance);
%             gama=abs(gama);
          gama= std_dev * randn(1, 1);
          if rand<0.5
                           OffspringDec(i+size(F1(1).ss,2),:)=FD(i,:)+DD*gama;
%              OffspringDec(i+size(F1(1).ss,2),:)=FN(selectedNumber,:)+DD*gama;
%               OffspringDec(i+size(F1(1).ss,2),:)=FN(selectedNumber,:)-DD*gama;
%             OffspringDec(i+size(F1(1).ss,2),:)=x_max-FD(i,:);
          else
%           OffspringDec(i+size(F1(1).ss,2),:)=FN(selectedNumber,:)+DD*gama;
                      OffspringDec(i+size(F1(1).ss,2),:)=FD(i,:)+rand*(x_max-FD(i,:));
%                       OffspringDec(i+size(F1(1).ss,2),:)=x_max-FD(i,:);
%                       OffspringDec(i+size(F1(1).ss,2),:)=FD(i,:)+DD*gama;
          end
      end 
%          OffspringDec11  =OffspringDec(1+size(F1(1).ss,2):size(F1(1).ss,2)+size(complementNumbers,2),:);
%              Lower = repmat(x_min,size(complementNumbers,2),1);
%         Upper = repmat(x_max,size(complementNumbers,2),1);
%         Site  = rand(size(complementNumbers,2),D) <pm;
%        mu    = rand(size(complementNumbers,2),D);
%        temp  = Site & mu<=0.5;
%      OffspringDec11(temp) =  OffspringDec11(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
%                          (1-( OffspringDec11(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(yita2+1)).^(1/(yita2+1))-1);
%     temp = Site & mu>0.5; 
%      OffspringDec11(temp) =  OffspringDec11(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
%                          (1-(Upper(temp)- OffspringDec11(temp))./(Upper(temp)-Lower(temp))).^(yita2+1)).^(1/(yita2+1)));
%        
%       OffspringDec(1+size(F1(1).ss,2):size(F1(1).ss,2)+size(complementNumbers,2),:)     =OffspringDec11;
     else
         if size(F1(1).ss,2)==1
          FN=ParentDec;
          SORS=0.00000000001;
          
     else
         for i=1:size(F1(1).ss,2)
         FN=ParentDec;
         distances = pdist2(FN(i,:), FN);
         SOR=sort(distances);
         SORS(i)=SOR(2);
         end
     end
         DAV=sum(SORS)/size(F1(1).ss,2);
         
         for i=1:size(F1(1).ss,2)
         if SORS(i)>DAV
           if size(F1(1).ss,2)==1
             if rand<0.5
             DD=(FN(i,:)-x_max)/norm(FN(i,:)-x_max);
             else
             DD=(FN(i,:)-x_min)/norm(FN(i,:)-x_min);
             end    
%                 gama=abs(gama);
          gama=  randn(1, 1);
          OffspringDec(i,:)=FN(i,:)+DD*gama;
           else
            NUM=1:size(F1(1).ss,2);
           HBNUM = setdiff(NUM, i);
           numElements = length(HBNUM);
          % 生成一个随机索引，范围从 1 到数字集合的大小
          randomIndex = randi(numElements);
          selectedNumber = HBNUM(randomIndex);
          DD=(FN(i,:)-FN(selectedNumber,:))/norm(FN(i,:)-FN(selectedNumber,:));
          for k=1:size(F1(1).ss,2)
             namda(k)=(FN(i,:)-FN(k,:))*DD';      
          end
         variance = var(namda);
         std_dev = sqrt(variance);
          
         gama= std_dev * randn(1, 1);
%            gama=abs(gama);
         OffspringDec(i,:)=FN(i,:)+DD*gama;
           end
         else
             if rand<0.5
             DD=(FN(i,:)-x_max)/norm(FN(i,:)-x_max);
             else
             DD=(FN(i,:)-x_min)/norm(FN(i,:)-x_min);
             end
          for k=1:size(F1(1).ss,2)
             namda(k)=(FN(i,:)-FN(k,:))*DD';      
          end
         variance = var(namda);
         std_dev = sqrt(variance);
            
         gama= std_dev * randn(1, 1);
%          gama=abs(gama);
         OffspringDec(i,:)=FN(i,:)+DD*gama;
         end

         end 
     end
% OffspringDec
    %% Polynomial mutation
% % % %     Lower = repmat(x_min,N,1);
% % % %     Upper = repmat(x_max,N,1);
% % % %     Site  = rand(N,D) <pm;
% % % %     mu    = rand(N,D);
% % % %     temp  = Site & mu<=0.5;
% % % %     OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
% % % %                          (1-(OffspringDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(yita2+1)).^(1/(yita2+1))-1);
% % % %     temp = Site & mu>0.5; 
% % % %     OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
% % % %                          (1-(Upper(temp)-OffspringDec(temp))./(Upper(temp)-Lower(temp))).^(yita2+1)).^(1/(yita2+1)));
% % % %        
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
    for i=1:N
        [OffspringDec(i,(x_num+1):(x_num+f_num)),~]=object_fun(OffspringDec(i,1:x_num),f_num,x_num,fun);
    end
% %模拟二进制交叉与多项式变异
% [pop,~]=size(chromo_parent);
% suoyin=1;
% for i=1:pop
%    %%%模拟二进制交叉
%    %初始化子代种群
%    %随机选取两个父代个体
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
% %    %确定两个父代个体不是同一个
% %    while isequal(chromo_parent(parent_1,:),chromo_parent(parent_2,:))
% %        parent_2=round(pop*rand(1));
% %        if(parent_2<1)
% %            parent_2=1;
% %        end
% %    end
%    chromo_parent_1=chromo_parent(parent_1,:);
%    chromo_parent_2=chromo_parent(parent_2,:);
%    off_1=chromo_parent_1;
%    off_2=chromo_parent_2;%这是不是就有问题了呢，好吧，我已经改成了2
%    if(rand(1)<pc)
%        %进行模拟二进制交叉
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
%            %使子代在定义域内
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
%        %计算子代个体的目标函数值
%        [off_1(1,(x_num+1):(x_num+f_num)),~]=object_fun(off_1,f_num,x_num,fun);
%        [off_2(1,(x_num+1):(x_num+f_num)),~]=object_fun(off_2,f_num,x_num,fun);
%    end
%    %%%多项式变异
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
%            %使子代在定义域内
%            if(off_1(j)>x_max(j))
%                off_1(j)=x_max(j);
%            elseif(off_1(j)<x_min(j))
%                off_1(j)=x_min(j);
%            end
%        end
%        %计算子代个体的目标函数值
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
%            %使子代在定义域内
%            if(off_2(j)>x_max(j))
%                off_2(j)=x_max(j);
%            elseif(off_2(j)<x_min(j))
%                off_2(j)=x_min(j);
%            end
%        end
%        %计算子代个体的目标函数值
%        [off_2(1,(x_num+1):(x_num+f_num)),~]=object_fun(off_2,f_num,x_num,fun);
%    end
%    off(suoyin,:)=off_1;
%    off(suoyin+1,:)=off_2;
%    suoyin=suoyin+2;
% end
% pick = randperm(2*pop, pop)'; % 从 2*NP 中随机选出 NP 个
%    chromo_offspring= off(pick, :);
 chromo_offspring=OffspringDec;
end
