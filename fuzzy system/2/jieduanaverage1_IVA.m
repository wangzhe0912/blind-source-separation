function [ newaveragey ,newaverage2y ] = jieduanaverage1_IVA( iter,Y,oldaveragey,oldaverage2y )
%计算y的平均值的更新公式
%lamuda是遗忘因子；
%k是当前时刻点；
%y是输出信号；
%averagey是y的平均值；
lamuda=0.999;
newaveragey=lamuda*(iter-1)/iter*oldaveragey+1/iter*Y(:,iter,:);
newaverage2y=lamuda*(iter-1)/iter*oldaverage2y+1/iter*(Y(:,iter,:).^3);

end

