function [ guanxi ] = xiangguanxing( y )
%计算相关性
normy=y/norm(y);
guanxi=0;
for i=1:4
    fzhi = Dhanshu( i,normy )+HDhanshu( i,normy );
    guanxi=guanxi+fzhi;
end

