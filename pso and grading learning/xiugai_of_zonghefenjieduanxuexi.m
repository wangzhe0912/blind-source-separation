 %线上处理，参数为6维
clear all
close all
clc
tic

for runtimes=1:2
Fs=10000;
T=1/Fs;
N=3000;

%录入源信号并生成混合信号
s = zeros(4,N);
for i=1:N
    s(1,i)=sin(2*pi*800*i*T)*sin(2*pi*25*i*T);    
    s(2,i)=sin(2*pi*300*i*T+6*cos(2*pi*60*i*T));
    s(3,i)=2*rand(1,1)-1;
    s(4,i)=sign(cos(2*pi*155*i*T)); 
end 
[B,N]=size(s);                                    %得到源信号的尺寸
figure                                            %源信号图
for i=1:B  
subplot(B,1,i);plot(s(i,:));
end
A=2*rand(B)-1          ;                              %混合矩阵,[-a,a]的均匀分布
x=A*s;                                            %观测信号
figure                                            %观察信号图
for i=1:B
    subplot(B,1,i);plot(x(i,:));
end


%初始化分离矩阵，粒子群相关参数及相关性参数
W=0.5*eye(B);


particle_number = 10;        % 粒子数
Dimension=1;                 % 粒子维数
q(:,1) = 0.01+0.03*rand(particle_number,1);        % 在[0.1,0.4]中随机产生的旋转角度作为粒子的初始位置
q(:,2) = 0*rand(particle_number,1); % 初始化粒子的移动速度 
q(:,3) = 0*rand(particle_number,1);                  % q为20*7
q_xsf = q;
gbest_xsf=0.04*rand;
pbest_xsf = 0.04*rand(particle_number,1);                        % 初始化个体最优位置 
vmax = 0.04;                % 最大移动速度0.05
num_sum =1;              % 迭代次数
gene=0.24;  % 惯性权重，权重递减  

oldaveragey=zeros(B,1);
oldaverage2y=zeros(B,1);
Cij=zeros(B,B);HCij=zeros(B,B);Cii=zeros(1,B);HCii=zeros(1,B);


jishu=0;   %用于记录用粒子群确定步长的次数
y(:,1)=x(:,1);
%开始逐个点分离信号
for iter=2:N
    %计算分离信号的相关性
    y(:,iter)=W *x(:,iter);
    [ newaveragey ,newaverage2y ] = jieduanaverage1( iter,y,oldaveragey,oldaverage2y );
    [ Cij , HCij ] = jieduancovyiyj( Cij , HCij ,  iter , y , newaveragey ,oldaveragey ,newaverage2y ,oldaverage2y);
    [ Cii , HCii ] = jieduancovyi( Cii , HCii,  iter , y , newaveragey ,oldaveragey ,newaverage2y ,oldaverage2y );
    oldaveragey=newaveragey;oldaverage2y=newaverage2y;
    lamuda =0.999;
    for i=1:4
         for j=1:4
             SC(i,j)=Cij(i,j)/(Cii(i)*Cii(j)).^.5;
             HC(i,j)=HCij(i,j)/(HCii(i)*Cii(j)).^.5;
         end
    end
    for i=1:4
         for j=1:4
              Dij(i,j)=max([abs(SC(i,j)),abs(HC(i,j)),abs(HC(j,i))]);
         end    
    end
    Dij=Dij-eye(4);
    for i=1:4
         Di(i)=max(Dij(i,:));
    end
    D=max(Di);


    %根据相关性分别选择粒子群或函数来确定步长
    if D>0.25
         %记录粒子群得到的数据
         jishu=jishu+1;
   
          kur_gbest_xsf = 0;                            % 初始化全局（群体）最优峰度值
          kur_pbest_xsf = zeros(particle_number,1);
          q(:,1) =  0.04*rand(particle_number,1);        % 在[0，0.2]中随机产生的旋转角度作为粒子的初始位置
          q_xsf = q;
                  
          for j = 1:particle_number            % 各粒子 n=20
               for kk = 1:Dimension        % 粒子各维 N=3
                    pre = q_xsf(j,kk);                     % 第j个粒子位置坐标中第k维的当前值
                    sub1 = pbest_xsf(j,kk) - pre;          % 此粒子（局部最优位置坐标中第k维值-当前位置的第k维的当前值）
                    sub2 = gbest_xsf(kk) - pre;            % 全局最优位置坐标中第k维值-此粒子当前位置的第k维的当前值
                    prev = q_xsf(j,kk+Dimension);                  % 第j个粒子，第k维的当前速度
                    tempv =sign(rand-0.5)* (gene * prev + 0.5 * rand(1) * sub1 + 0.5 * rand(1) * sub2);   % 更新粒子的速度
                                      % 粒子每一维的运动速度V都被限制在[-Vmax,Vmax]之间
                    if tempv > vmax
                          q_xsf(j,kk + Dimension) = vmax;
                    elseif tempv < -1 * vmax
                          q_xsf(j,kk + Dimension) = -vmax;
                    else
                          q_xsf(j,kk + Dimension) = tempv;
                    end
                end
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % 粒子位置的更新
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 for kk = 1:Dimension
                       xx = q_xsf(j,kk) + q_xsf(j,kk + Dimension);        % Xi(n+1)=Xi(n)+Vi(n)
                       if xx >0.04
                            xx = 0.04;
                       end
                       if xx < 0.01 
                            xx = 0.01;
                       end
                            q_xsf(j,kk) = xx;
                  end                   
                       buchang=diag(q_xsf(j,1:Dimension)); 
                       kur_y1 = zonghanshu( buchang,x,iter,W ,Cij,HCij,Cii,HCii,oldaveragey,oldaverage2y,Di) ;                 % 计算分离出的信号的峰度（根据公式6）
                       kur_pbest_xsf(j) =zonghanshu(diag(pbest_xsf(j,:)) ,x,iter,W ,Cij,HCij,Cii,HCii,oldaveragey,oldaverage2y,Di) ;   
                       kur_gbest_xsf=zonghanshu(diag(gbest_xsf) ,x,iter,W ,Cij,HCij,Cii,HCii,oldaveragey,oldaverage2y,Di) ;   
                       q_xsf(j,Dimension*2+1) = kur_y1;                      % 保存粒子当前峰度值
                  if kur_y1 <  kur_pbest_xsf(j)              % 搜索个体最优位置
                        kur_pbest_xsf(j) = kur_y1;
                        pbest_xsf(j,:) = q_xsf(j,1:Dimension);
                  end
                  if kur_y1 <  kur_gbest_xsf                 % 搜索全局（群体）最优位置
                        kur_gbest_xsf = kur_y1;               % 全局最优峰度值
                        gbest_xsf = q_xsf(j,1:Dimension);             % 全局最优位置
                  end
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          end
     
          buchang=diag(gbest_xsf);
          buchangjilu(jishu,:)=gbest_xsf;

    else
         buchang=zeros(4,4);
         for ii=1:4
              if Di(i)>0.2
                   buchang(ii,ii)=0.016;
              elseif Di(i)>0.05
                   buchang(ii,ii)=0.00275+0.05*(Di(ii)-0.05)^0.7;
              else
                   buchang(ii,ii)=0.135*Di(ii)^1.3;
              end
          end
    end
 
     W = W + buchang*( eye(B) - ((W *x(:,iter)).^2).*sign(W *x(:,iter))*3*tanh(10*W *x(:,iter))' )*W;
     
     E(iter)=crosstalking ( W,A );
     
end
y=real(y);

figure                                           %估计源信号图
for nn=1:B
    subplot(B,1,nn);plot(y(nn,:));
end

%截取第3001-4000的点用于衡量分离性能
    ss=s(:,2001:3000);
    yy=y(:,2001:3000);
for i=1:4
    for j=1:4
        a=corrcoef(ss(i,:),yy(j,:));   %计算源信号与分离信号的相关系数
        xingneng(i,j)=a(1,2);
    end
end
        abs(xingneng)
        
        
%调整分离信号的顺序与正负值
[B,shunxu]=max(abs(xingneng));   %找出分离信号与源信号的对应关系
for i=1:4
    yyy(shunxu(i),:)=yy(i,:);    %重新分配，使得源信号与分离信号顺序上对应
end
for i=1:4                        %判断源信号与分离信号的正负关系并修正
    ans=corrcoef(yyy(i,:),ss(i,:));    
    if ans(1,2)<0
        yyy(i,:)=-yyy(i,:);
    end
end
    figure                       %绘制重新调整后的分离信号，并取其3001:4000的点
    for i=1:4
        subplot(4,1,i)
        plot(yyy(i,:));
    end

%%计算信噪比
for i=1:4                         %对源信号与分离信号进行归一化
    ss(i,:)=ss(i,:)/norm(ss(i,:));
    yyy(i,:)=yyy(i,:)/norm(yyy(i,:));
end
[B,changdu]=size(ss);                    %得到源信号的尺寸
sums=zeros(1,4);                   %赋初值
sume=zeros(1,4);
for i=1:B                          %计算信噪比
    ss2(i,:)=ss(i,:).^2;
    error1(i,:)=(yyy(i,:)-ss(i,:)).^2;
    for j=1:changdu
        sums(i)=sums(i)+ss2(i,j);
    end
    for j=1:changdu
        sume(i)=sume(i)+error1(i,j);
    end
    snr(i)=10*log10(sums(i)/sume(i));
end
snr                                %显示信噪比  
figure;
plot(E)
Etongji(runtimes,:)=E;
chouyangE=E(2901:3000);
meanE=mean(chouyangE)
meanEtongji(runtimes)=meanE;
stdE=std(chouyangE)
stdEtongji(runtimes)=stdE;
end
toc