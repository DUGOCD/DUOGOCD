function chromo = initialize1( pop,f_num,x_num,x_min,x_max )
%   种群初始化
for i=1:pop
   for j=1:x_num
       chromo(i,j)=x_min(j)+(x_max(j)-x_min(j))*rand(1);
   end
   chromo(i,x_num+1)=i;
end
end