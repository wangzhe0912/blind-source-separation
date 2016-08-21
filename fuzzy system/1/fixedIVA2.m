clc;
clear all;
close all;
tic
%% 录入源信号并生成混合信号
%读取信号
load('70s_signals.mat');
begin=fix(rand*56000);

s1=sss(2,begin+1:begin+114000);
s2=sss(4,begin+1:begin+114000);
s=[s1;s2];
%绘图
figure(1)
subplot(2,1,1)
plot(s1)
subplot(2,1,2)
plot(s2)
xlabel('采样点数');

%获取混合系数
RR=[1,2];
SS=[3,4];
hh=RIR2(RR,SS);
reveber_filter=1600;
h11=reshape(hh(1,1,:),1,reveber_filter);
h12=reshape(hh(1,2,:),1,reveber_filter);
h21=reshape(hh(2,1,:),1,reveber_filter);
h22=reshape(hh(2,2,:),1,reveber_filter);

%得到混合信号
x1=filter(h11,1,s1)+filter(h12,1,s2);
x2=filter(h21,1,s1)+filter(h22,1,s2);

x(1,:)=x1;
x(2,:)=x2;
figure(2)
subplot(2,1,1)
plot(x1)
subplot(2,1,2)
plot(x2)
xlabel('采样点数');


%% 对混合信号进行时频转换，得到X
[nsou, nn] = size(x);
if ~exist('nfft','var')|isempty(nfft), nfft = 256; end %采样点数256
%if ~exist('nfft','var')|isempty(nfft), nfft = 1024; end %采样点数1024
win = 2*hanning(nfft,'periodic')/nfft;
nol = fix(3*nfft/4);

for l=1:nsou,
    X(l,:,:) = conj(stft(x(l,:)', nfft, win, nol)');
end
clear x;
clear S;

%% 初始化分离矩阵，粒子群相关参数及相关性参数
N = size(X,2);
nfreq = size(X,3);
epsi = 1e-6;
pObj = Inf;

for k=1:nfreq,
    Wp(:,:,k) = eye(nsou);  
end

oldaveragey=zeros(nsou,1,nfreq);
oldaverage2y=zeros(nsou,1,nfreq);
Cij=zeros(nsou,nsou,nfreq);HCij=zeros(nsou,nsou,nfreq);Cii=zeros(1,nsou,nfreq);HCii=zeros(1,nsou,nfreq);

jishu=0;
etajilu=zeros(N,nfreq);
for k=1:nfreq,
    Y(:,1,k) = Wp(:,:,k)*X(:,1,k);    %输入一节数据点，以目前的分离矩阵处理当前一节数据
end
SIR1=zeros(15,2,2);
% Start iterative learning algorithm
for iter=2:N,
    
    for k=1:nfreq,
        Y(:,iter,k) = Wp(:,:,k)*X(:,iter,k);    %输入一节数据点，以目前的分离矩阵处理当前一节数据
    end
    
    [ newaveragey ,newaverage2y ] = jieduanaverage1_IVA( iter,abs(Y),oldaveragey,oldaverage2y );
    [ Cij , HCij ] = jieduancovyiyj_IVA( Cij , HCij ,  iter , abs(Y) , newaveragey ,oldaveragey ,newaverage2y ,oldaverage2y);
    [ Cii , HCii ] = jieduancovyi_IVA( Cii , HCii,  iter , abs(Y) , newaveragey , oldaveragey , newaverage2y , oldaverage2y );
    oldaveragey=newaveragey;
    oldaverage2y=newaverage2y;
    lamuda =0.999;
    for k=1:nfreq
        for i=1:nsou
             for j=1:nsou             
                 SC(i,j,k)=Cij(i,j,k)/(Cii(:,i,k).*Cii(:,j,k)).^.5;
                 HC(i,j,k)=HCij(i,j,k)/(HCii(:,i,k).*Cii(:,j,k)).^.5;
             end
         end
    end
    for k=1:nfreq
         for i=1:nsou
             for j=1:nsou             
                 Dijk(i,j,k)=max([abs(SC(i,j,k)),abs(HC(i,j,k)),abs(HC(j,i,k))]);
             end
         end    
    end
    for k=1:nfreq
         Dijk(:,:,k)=Dijk(:,:,k)-eye(nsou,nsou);
    end
    
    for k=1:nfreq
         for i=1:nsou
              Dik(i,k)=max(Dijk(i,:,k));
         end
    end
    D=max(Dik);   
    meanD(iter)=mean(D);
    
    Ssq = sum(abs(Y(:,iter,:)).^2,3).^.5;        %出现差异
    Ssq1 = (Ssq+epsi).^-1;

    for k=1:nfreq,
        % Calculate multivariate score function and gradients
        Phi = Ssq1.*Y(:,iter,k);
        duijiao=eye(nsou);
        R=Phi*Y(:,iter,k)';
        for i=1:nsou
            duijiao(i,i)=R(i,i);
        end
        dWp(:,:,k) = (duijiao - Phi*Y(:,iter,k)')*Wp(:,:,k);%对应于式（15）（16）

    end
    
    % Update unmixing matrices
    if iter==2
        xishu=ones(k,1);
    else
        xishu=0.5*xishu+0.5*0.5*squeeze(sum((abs(X(:,iter,:))).^2,1));
    end
    buchongxiang=xishu.^(-0.5);
    
    eta=0.3;
    %% 信号分离
    
    for k=1:nfreq
        Wp(:,:,k) = Wp(:,:,k) + eta*buchongxiang(k).*dWp(:,:,k);
    end
    %%最终处理
    if(mod(iter,125)==0)
        % Correct scaling of unmixing filter coefficients
        for k=1:nfreq,
              W(:,:,k) = Wp(:,:,k);
              W(:,:,k) = diag(diag(pinv(W(:,:,k))))*W(:,:,k);
        end

        % Calculate outputs
        for k=1:nfreq,
           YY(:,:,k) = W(:,:,k)*X(:,:,k);
        end

        % Re-synthesize the obtained source signals
        for k=1:nsou,
           y(k,:) = istft(conj(squeeze(YY(k,:,:))'), nn, win, nol)';
        end

        for i=1:nsou
            y(i,:)=y(i,:)/norm(y(i,:));
        end

        figure
        for i=1:2
             subplot(2,1,i)
             plot(y(i,:))
        end
        p=iter/125;
        %SIR1(p,:,:)=SIR(y(:,1:56000),s(:,1:56000),reveber_filter); %计算SIR（过于耗时）
           
    end
end

figure
plot(meanD)
toc