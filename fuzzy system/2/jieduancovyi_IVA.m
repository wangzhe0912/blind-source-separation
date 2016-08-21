function [ Cii , HCii ] = jieduancovyi_IVA( Cii , HCii,  iter , Y , newaveragey , oldaveragey , newaverage2y , oldaverage2y )
%¼ÆËãcov(yi)
B=size(Y,1);
nfreq=size(Y,3);
lamuda =0.999;
for k=1:nfreq
    for i=1:B
         Cii(:,i,k)=lamuda*(iter-1)/iter*(Cii(:,i,k)+(newaveragey(i,:,k)-oldaveragey(i,:,k)).^2)+1/iter*(Y(i,iter,k)-newaveragey(i,:,k)).^2;
         HCii(:,i,k)=lamuda*(iter-1)/iter*(HCii(:,i,k)+(newaverage2y(i,:,k)-oldaverage2y(i,:,k)).^2)+1/iter*((Y(i,iter,k).^3)-newaverage2y(i,:,k)).^2;
    end
end
end

