function [sumkurt] = sumkurt( y )
%求yi的俏度的绝对值之和
[n,k]=size(y);
sum=0;
for i=1:n
    sum=sum+abs(kurt(white(y(i,:))));
end
sumkurt=sum;
end

