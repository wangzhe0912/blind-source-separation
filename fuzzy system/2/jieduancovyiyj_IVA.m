function [ Cij , HCij ] = jieduancovyiyj_IVA( Cij , HCij ,  iter , Y , newaveragey ,oldaveragey ,newaverage2y ,oldaverage2y)
%¼ÆËãcov(yi,yj)
B=size(Y,1);
nfreq=size(Y,3);
lamuda =0.999;
for k=1:nfreq
    for i=1:B
        for j=1:B
          Cij(i,j,k)=lamuda*(iter-1)/iter*(Cij(i,j,k)+(newaveragey(i,:,k)-oldaveragey(i,:,k)).*(newaveragey(j,:,k)-oldaveragey(j,:,k)))+1/iter*(Y(i,iter,k)-newaveragey(i,:,k))*(Y(j,iter,k)-newaveragey(j,:,k));
          HCij(i,j,k)=lamuda*(iter-1)/iter*(HCij(i,j,k)+(newaverage2y(i,:,k)-oldaverage2y(i,:,k)).*(newaveragey(j,:,k)-oldaveragey(j,:,k)))+1/iter*((Y(i,iter,k).^3)-newaverage2y(i,:,k))*(Y(j,iter,k)-newaveragey(j,:,k));
        end
    end
end
end

