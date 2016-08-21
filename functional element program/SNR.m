function  [snr1,snr2]=SNR(s,y)
[B,P]=size(s);
sums=zeros(1,B);
sumy=zeros(1,B);
global shunxu;
for i=1:B
    s2(i,:)=s(i,:).^2;    
    y2(i,:)=y(i,:).^2; 
    for j=1:P
        sums(i)=sums(i)+s2(i,j);
        sumy(i)=sumy(i)+y2(i,j);
    end
end
if shunxu==1;
    error1(1,:)=(y(1,:)-s(1,:)).^2;
    error1(2,:)=(y(2,:)-s(2,:)).^2;
    sume1=0;
    sume2=0;
    for j=1:P
        sume1=sume1+error1(1,j);
        sume2=sume2+error1(2,j);
    end
    snr1=10*log10(sums(1)/sume1);
    snr2=10*log10(sums(2)/sume2);
end
if shunxu==0;
    error1(1,:)=(y(2,:)-s(1,:)).^2;
    error1(2,:)=(y(1,:)-s(2,:)).^2;
    sume1=0;
    sume2=0;
    for j=1:P
        sume1=sume1+error1(1,j);
        sume2=sume2+error1(2,j);
    end
    snr1=10*log10(sums(1)/sume1);
    snr2=10*log10(sums(2)/sume2);
    end
end
