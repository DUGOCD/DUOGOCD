function chromo = initialize2( complementNumbers,x_num,x_min,x_max )
%   种群初始化
for i=1:size(complementNumbers,2)
   for j=1:x_num
       chromo(i,j)=x_min(j)+(x_max(j)-x_min(j))*rand(1);
   end
   chromo(i,x_num+1)=complementNumbers(i);
end
end