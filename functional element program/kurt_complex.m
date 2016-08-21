function [ Kurt ] = kurt_complex( y )
%计算复数信号的峭度
y=white(y);
Kurt=mean(abs(y).^4)-2;

end

