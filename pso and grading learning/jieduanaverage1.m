function [ newaveragey ,newaverage2y ] = jieduanaverage1( k,y,oldaveragey,oldaverage2y )
%计算y的平均值的更新公式
%lamuda是遗忘因子；
%k是当前时刻点；
%y是输出信号；
%averagey是y的平均值；
lamuda=0.999;
newaveragey=lamuda*(k-1)/k*oldaveragey+1/k*y(:,k);
newaverage2y=lamuda*(k-1)/k*oldaverage2y+1/k*(y(:,k).^3);

end

