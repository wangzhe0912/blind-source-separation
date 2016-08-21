function [y] = shunshipso(x,s)
%  初始化粒子群
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 20;        % 粒子数
N=1;
%lieshu = N * N * 2 + 1;                       % 19 N=3是信号的个数
%p = rand(n,lieshu);                          % p为产生w的随机矩阵与初始速度的组合 ; 
                                              % Uniformly distributed(均匀分布）
% N位为粒子的初始位置，中间N位为初始速度，后面为适值(峰度值)
q(:,1) =  2*pi*(rand(n,1) - 0.5);        % 在[-pi，pi]中随机产生的旋转角度作为粒子的初始位置
q(:,2) = 0.001 * (rand(n,1) - 0.5); % 初始化粒子的移动速度 
q(:,3) = zeros(n,1);                  % q为20*7
q_xsf = q;

% 初始化pbest_xsf和gbest_xsf及其适应度
pbest_xsf = rand(n,1);                        % 初始化个体最优位置
kur_gbest_xsf = 0;                            % 初始化全局（群体）最优峰度值
gbest_xsf=zeros(1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 初始化，确定每一个粒子的初始位置和峰度
% 找到初始状态下的全局最优粒子和最大峰度值      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n     %n=20粒子数
    ww=cir2(q_xsf(i,1)); % 确定分离矩阵，q(i,1),q(i,2),q(i,3)为粒子，随机产生
    pbest_xsf(i,:) = q_xsf(i,1);            % 存放粒子的位置
    y1 = ww * x;   
    kur_y1 = kurtosis(y1);                    % kurtosis 求峭度
    %kur_y1 = mi(y1,ww,x); 
    
    q_xsf(i,2) = kur_y1;                      % 存放第i个粒子的峰度值
    kur_pbest_xsf(i) = kur_y1;
    if kur_pbest_xsf(i) > kur_gbest_xsf
        gbest_xsf = pbest_xsf(i,:);
        kur_gbest_xsf = kur_pbest_xsf(i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 当为0.09时输出为465.6466，0.03时为542.4920
gene_max = 0.8;            % 最大权重
gene_min = 0.3;            % 最小权重
vmax = 0.7;                % 最大移动速度0.05
num_sum =50;              % 迭代次数
step_max = 0.5;            % 最大加权因子 （用于整体速度加权）
step_min = 0.5;            % 最小加权因子
for i = 1:num_sum          % 迭代步数
    gene = gene_max - (gene_max - gene_min) * i / num_sum;  % 惯性权重，权重递减
    step = step_max - (step_max - step_min) * i / num_sum;       
    for j = 1:n            % 各粒子 n=20
        ww=cir2(q_xsf(j,1)); % 计算第j个粒子在目前位置所代表的解混矩阵w
        y1 = ww * x;                              % 得到由此w计算出的分离结果       
        kur_y1 = kurtosis(y1);                    % 计算分离出的信号的峰度（根据公式6）
        %kur_y1 = mi(y1,ww,x);
        
        q_xsf(j,3) = kur_y1;                      % 保存粒子当前峰度值
        if kur_y1 > kur_pbest_xsf(j)              % 搜索个体最优位置
            kur_pbest_xsf(j) = kur_y1;
            pbest_xsf(j,:) = q_xsf(j,1);
        end
        if kur_y1 > kur_gbest_xsf                 % 搜索全局（群体）最优位置
            kur_gbest_xsf = kur_y1;               % 全局最优峰度值
            gbest_xsf = q_xsf(j,1);             % 全局最优位置
        end
        k_best_xsf(i) = kur_gbest_xsf;            % 记录每次迭代的全局最优峰度值
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 粒子速度的更新
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1        % 粒子各维 N=3
            % x_bc(i,k) = std(q_xsf(:,k))^3;      % std:标准偏差  ？？？？？？
            pre = q_xsf(j,k);                     % 第j个粒子位置坐标中第k维的当前值
            sub1 = pbest_xsf(j,k) - pre;          % 此粒子（局部最优位置坐标中第k维值-当前位置的第k维的当前值）
            sub2 = gbest_xsf(k) - pre;            % 全局最优位置坐标中第k维值-此粒子当前位置的第k维的当前值
            prev = q_xsf(j,k+1);                  % 第j个粒子，第k维的当前速度
            tempv = step * (gene * prev + 2 * rand(1) * sub1 + 2 * rand(1) * sub2);   % 更新粒子的速度
            g_bc(i,j,k) = tempv;                  % 第i次迭代中，第j个粒子的第k维新速度
            % 粒子每一维的运动速度V都被限制在[-Vmax,Vmax]之间
            if tempv > vmax
                q_xsf(j,k + 1) = vmax;
            elseif tempv < -1 * vmax
                q_xsf(j,k + 1) = -vmax;
            else
                q_xsf(j,k + 1) = tempv;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 粒子位置的更新
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:1
             xx = q_xsf(j,k) + q_xsf(j,k + 1);        % Xi(n+1)=Xi(n)+Vi(n)
             if xx > 2*pi
                 xx = xx - 2 * pi;
             end
             if xx < -2 * pi
                 xx = xx + 2 * pi;
             end
            q_xsf(j,k) = xx;
        end 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    q_pace(i,:,:) = q_xsf;
    qq_gbest_xsf(i,:) = gbest_xsf;
    kur_gbest_xsf;
end

% 信号分离
w=cir2(gbest_xsf)   % 找到的最优粒子位置，计算出分离矩阵
y = w * x;                                        % 得到最终分离信号 

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
%qq=[s(3,:)' y(3,:)'];
fprintf('信源与恢复信号相关性');
format long
cc=corrcoef(bb)  % 源信号与恢复的信号的相关系数
ee=corrcoef(dd)
cc2=corrcoef(bb1)  % 源信号与恢复的信号的相关系数
ee2=corrcoef(dd1)
%ff=corrcoef(qq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global shunxu;
if abs(cc(1,2))>abs(cc2(1,2));
    shunxu=1
else shunxu=0
end

figure(5)
hold on
plot(k_best_xsf,'g');
xlabel('迭代步数');
ylabel('适度值');

hold off