function [ E ] = crosstalking ( W,A )
%jisuanCROSSTALING
P=W*A;
[num,b]=size(A);
sum1=0;sum2=0;
for i=1:num
    for j=1:num
       sum1=sum1+ abs(P(i,j))/max(abs(P(i,:))); sum2=sum2+ abs(P(j,i))/max(abs(P(:,i))); 
    end
    sum1=sum1-1;sum2=sum2-1;
end

E=(sum1+sum2)/num;
end

