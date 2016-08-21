clc;
clear all;
close all;
tic
%% 录入源信号并生成混合信号
%读取信号
load('70s_signals.mat');
begin=fix(rand*56000);
begin=0;
s1=sss(1,begin+1:begin+224000);
s2=sss(2,begin+1:begin+224000);
s3=sss(3,begin+1:begin+224000);
s4=sss(4,begin+1:begin+224000);
s=[s1;s2;s3;s4];
%绘图
size1=size(s,1);
length=size(s,2);
figure(1)
for i=1:size1
    subplot(size1,1,i)
    plot(s(i,:))
end
xlabel('采样点数');

%获取混合系数
RR=[1,2,3,4];
SS=[1,2,3,4];
hh=RIR2(RR,SS);
reveber_filter=1600;

x=zeros(size1,length);
for i=1:size1
    for j=1:size1
    x(i,:)=x(i,:)+filter(reshape(hh(i,j,:),1,reveber_filter),1,s(j,:));
    end
end

figure(2)
for i=1:size1
    subplot(size1,1,i)
    plot(x(i,:))
end
xlabel('采样点数');


%% 对混合信号进行时频转换，得到X
[nsou, nn] = size(x);
%if ~exist('nfft','var')|isempty(nfft), nfft = 1024; end
if ~exist('nfft','var')|isempty(nfft), nfft = 256; end
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

for k=1:nfreq,
    Y(:,1,k) = Wp(:,:,k)*X(:,1,k);    %输入一节数据点，以目前的分离矩阵处理当前一节数据
end

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

    
   %% 信号分离

    eta=zeros(1,nfreq);
    %% 基于模糊推断系统的步长选择
    fuzzy_function = readfis('fuzzy_function1');
    for k=1:nfreq
       eta(k) = evalfis([meanD(iter),D(k)],fuzzy_function);
    end
    eta=1*eta;
    for k=1:nfreq
        Wp(:,:,k) = Wp(:,:,k) + eta(k)*buchongxiang(k).*dWp(:,:,k);
    end


end

% Correct scaling of unmixing filter coefficients
for k=1:nfreq,
    W(:,:,k) = Wp(:,:,k);
    W(:,:,k) = diag(diag(pinv(W(:,:,k))))*W(:,:,k);
end

% Calculate outputs
for k=1:nfreq,
    finalY(:,:,k) = W(:,:,k)*X(:,:,k);
end

% Re-synthesize the obtained source signals
for k=1:nsou,
    finaly(k,:) = istft(conj(squeeze(finalY(k,:,:))'), nn, win, nol)';
    y(k,:) = istft(conj(squeeze(Y(k,:,:))'), nn, win, nol)';
end

for i=1:nsou
    finaly(i,:)=finaly(i,:)/norm(finaly(i,:));
    y(i,:)=y(i,:)/norm(y(i,:));
end

figure(3)
for i=1:size1
    subplot(size1,1,i)
    plot(finaly(i,:))
end
xlabel('采样点数');

figure(4)
for i=1:size1
    subplot(size1,1,i)
    plot(y(i,:))
end
xlabel('采样点数');

%SIR1=SIR(finaly(:,1:56000),s(:,1:56000),reveber_filter);
figure
plot(meanD)
toc