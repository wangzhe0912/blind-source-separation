function [y] = newjuanjipso(x,s)
n = 20;         % 粒子数
N = 40;         % 粒子维数
L=30;           % 互累积量阶数
qmax = 1;qmin = -1;%粒子位置范围
for i = 1:n
    for j=1:N
    q_xsf(i,j)=(qmax-qmin)*rand+qmin;%粒子位置初始化
    end
end
for i=1:n
   for j=1:N
        v(i,j)=0.001 * (rand-0.5) ;  % 粒子速度随机初始化
   end
   
end
pbest_xsf = rand(n,N);                        % 初始化粒子个体最优位置
gbest_xsf = rand(1,N);
kur_gbest_xsf = 0;                            % 初始化全局（群体）最优峰度值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 初始化，确定每一个粒子的初始位置和峰度
% 找到初始状态下的全局最优粒子和最大峰度值      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n
    ww=reshape(q_xsf(i,:),4,10);       % 由粒子位置确定分离矩阵
    pbest_xsf(i,:) = q_xsf(i,1:40);     % 存放粒子的位置
   y(1,:) = filter(ww(1,:),1,x(1,:))+filter(ww(2,:),1,x(2,:));
   y(2,:) = filter(ww(3,:),1,x(1,:))+filter(ww(4,:),1,x(2,:));
   y = white(y);                    % 白化——很关键
   y1=y(1,:);
   y2=y(2,:);
  %  kur_y1 = leijiliang_5(x1,x2,ww,L); 
   kur_y1=leijiliang(y1,y2,L);
    %kur_y1 = leijiliang_yi(x1,x2,ww,y1,y2,L);% 求四阶互累积量
    p(i) = kur_y1;                     % 存放第i个粒子的峰度值
    kur_pbest_xsf(i) = kur_y1;         % 把峰度值给个体最优峰度值
    if kur_pbest_xsf(i) > kur_gbest_xsf
        gbest_xsf = pbest_xsf(i,:);    %找到群最优位置
        kur_gbest_xsf = kur_pbest_xsf(i);%群最优峰度值
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSO优化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gene_max = 0.8;            % 最大权重
gene_min = 0.3;            % 最小权重
vmax = 0.7;                % 最大移动速度0.05
num_sum =60;               % 迭代次数
step_max = 0.5;            % 最大加权因子 （用于整体速度加权）
step_min = 0.5;            % 最小加权因子
for i = 1:num_sum          % 迭代步数
    gene = gene_max - (gene_max - gene_min) * i / num_sum;  % 惯性权重，权重递减
    step = step_max - (step_max - step_min) * i / num_sum;  %基本没用    
    for j = 1:n                      % 各粒子
      ww=reshape(q_xsf(j,:),4,10); % 计算第j个粒子在目前位置所代表的解混矩阵w
      y(1,:) = filter(ww(1,:),1,x(1,:))+filter(ww(2,:),1,x(2,:));
      y(2,:) = filter(ww(3,:),1,x(1,:))+filter(ww(4,:),1,x(2,:));
      y = white(y);                    % 白化——很关键
      y1=y(1,:);
      y2=y(2,:);
     kur_y1 = leijiliang(y1,y2,L);             % 求四阶互累积量
       % kur_y1 = leijiliang_5(x1,x2,ww,L); 
      %   kur_y1=abs(cum31(y1,y2)+cum22(y1,y2));
     % kur_y1 = leijiliang_yi(x1,x2,ww,y1,y2,L);
        p(j) = kur_y1;                      % 保存粒子当前峰度值
        if kur_y1 > kur_pbest_xsf(j)              % 搜索个体最优位置
            kur_pbest_xsf(j) = kur_y1;
            pbest_xsf(j,:) = q_xsf(j,1:N);        % 个体最优位置
        end
        if kur_y1 > kur_gbest_xsf                 % 搜索全局（群体）最优位置
            kur_gbest_xsf = kur_y1;               % 全局最优峰度值
            gbest_xsf = q_xsf(j,1:N);             % 全局最优位置
        end
        k_best_xsf(i) = kur_gbest_xsf;            % 记录每次迭代的全局最优峰度值
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 粒子速度的更新
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:N        % 粒子各维
            pre = q_xsf(j,k);                     % 第j个粒子位置坐标中第k维的当前值
            sub1 = pbest_xsf(j,k) - pre;          % 此粒子（局部最优位置坐标中第k维值-当前位置的第k维的当前值）
            sub2 = gbest_xsf(k) - pre;            % 全局最优位置坐标中第k维值-此粒子当前位置的第k维的当前值
            prev = v(j,k);                        % 第j个粒子，第k维的当前速度
            tempv = step * (gene * prev + 2 * rand(1) * sub1 + 2 * rand(1) * sub2);   % 更新粒子的速度
            g_bc(i,j,k) = tempv;                  % 第i次迭代中，第j个粒子的第k维新速度（没啥用）
            % 粒子每一维的运动速度V都被限制在[-Vmax,Vmax]之间
            if tempv > vmax
                v(j,k ) = vmax;
            elseif tempv < -1 * vmax
                v(j,k ) = -vmax;
            else
                v(j,k ) = tempv;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 粒子位置的更新
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for k = 1:N
         q_xsf(j,k)= q_xsf(j,k) + v(j,k );
           if q_xsf(j,k) > qmax
               q_xsf(j,k) = qmax;
           end
           if q_xsf(j,k) < qmin
               q_xsf(j,k) = qmin;
           end
      end
    end
    q_pace(i,:,:) = q_xsf;              %每一次迭代最后一个粒子的位置
    qq_gbest_xsf(i,:) = gbest_xsf;      %群最优位置
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % 信号分离   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=reshape(gbest_xsf,4,10);             % 根据找到的最优粒子位置，计算出分离矩阵
 y(1,:) = filter(w(1,:),1,x(1,:))+filter(w(2,:),1,x(2,:));
 y(2,:) = filter(w(3,:),1,x(1,:))+filter(w(4,:),1,x(2,:));
  y = white(y);     % 得到最终分离信号
                     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 绘制最后的输出图    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hdl=figure('Name','(3)分离后的信号','NumberTitle','off','MenuBar','none','Position',[1000 300 400 400]);
for i = 1:2
   subplot(2,1,i)
   plot(y(i,:));
   xlabel('采样点');
   ylabel('幅度');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 分离性能分析
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
gplotmatrix(s',y')   %误差，即w*A  
bb=[s(1,:)' y(1,:)'];
dd=[s(2,:)' y(2,:)'];
bb1=[s(1,:)' y(2,:)'];
dd1=[s(2,:)' y(1,:)'];

fprintf('信源与恢复信号相关性');
format long
cc=corrcoef(bb)  % 源信号与恢复的信号的相关系数
ee=corrcoef(dd)
cc2=corrcoef(bb1)  % 源信号与恢复的信号的相关系数
ee2=corrcoef(dd1)
%ff=corrcoef(qq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
hold on
plot(k_best_xsf,'g');
xlabel('迭代步数');
ylabel('适度值');

hold off