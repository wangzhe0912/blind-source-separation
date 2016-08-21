function [ f ] = ZONGHANSHU( buchang,x,k,W ,Cij,HCij,Cii,HCii,oldaveragey,oldaverage2y,Di)
%整体函数，用于体现W与f的关系
     B=4;
     W = W + buchang*( eye(B) - ((W *x(:,k)).^2).*sign(W *x(:,k))*3*tanh(10*W *x(:,k))' )*W;
     y(:,k)=W *x(:,k);
     [ newaveragey ,newaverage2y ] = jieduanaverage1( k,y,oldaveragey,oldaverage2y );
     [ Cij , HCij ] = jieduancovyiyj( Cij , HCij ,  k , y , newaveragey ,oldaveragey ,newaverage2y ,oldaverage2y);
     [ Cii , HCii ] = jieduancovyi( Cii , HCii,  k , y , newaveragey , oldaveragey , newaverage2y , oldaverage2y );
     for i=1:4
           for j=1:4
                SC(i,j)=Cij(i,j)/(Cii(i)*Cii(j)).^.5;
                HC(i,j)=HCij(i,j)/(HCii(i)*Cii(j)).^.5;
           end
     end
     for i=1:4
           for j=1:4
                Dij(i,j)=max([abs(SC(i,j)),abs(HC(i,j)),abs(HC(j,i))]);
           end    
     end
     Dij=Dij-eye(4);
     for i=1:4
           newDi(i)=max(Dij(i,:));
     end
     f=objection(newDi,Di);
end

