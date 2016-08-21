function [ Cii , HCii ] = jieduancovyi( Cii , HCii,  k , y , newaveragey , oldaveragey , newaverage2y , oldaverage2y )
%¼ÆËãcov(yi)
B=4;
lamuda =0.999;
for i=1:B
     Cii(i)=lamuda*(k-1)/k*(Cii(i)+(newaveragey(i)-oldaveragey(i)).^2)+1/k*(y(i,k)-newaveragey(i)).^2;
     HCii(i)=lamuda*(k-1)/k*(HCii(i)+(newaverage2y(i)-oldaverage2y(i)).^2)+1/k*((y(i,k).^3)-newaverage2y(i)).^2;
end
end

