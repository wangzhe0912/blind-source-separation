function [ Cij , HCij ] = jieduancovyiyj( Cij , HCij ,  k , Y , newaveragey ,oldaveragey ,newaverage2y ,oldaverage2y)
%¼ÆËãcov(yi,yj)
B=size(Y,1);
lamuda =0.999;
for i=1:B
    for j=1:B
          Cij(i,j)=lamuda*(k-1)/k*(Cij(i,j)+(newaveragey(i)-oldaveragey(i))*(newaveragey(j)-oldaveragey(j)))+1/k*(Y(i,k)-newaveragey(i))*(Y(j,k)-newaveragey(j));
          HCij(i,j)=lamuda*(k-1)/k*(HCij(i,j)+(newaverage2y(i)-oldaverage2y(i))*(newaveragey(j)-oldaveragey(j)))+1/k*((Y(i,k).^3)-newaverage2y(i))*(Y(j,k)-newaveragey(j));
    end
end
end

