function [ HD ] = HDhanshu( i,y )
%求高阶相关
[a,b]=size(y);
HD=0;
for j=1:a
A=corrcoef(y(i,:).^2+y(i,:).^3,y(j,:).^2+y(j,:).^3);
HD=HD+abs(A(1,2)^2);  
end
HD=HD-1;
end
