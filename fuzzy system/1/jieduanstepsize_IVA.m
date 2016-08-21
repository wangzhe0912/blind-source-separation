function [ buchang, Dij, Di, D ,SC ,HC] = jieduanstepsize(  k , R , P , S , Q )
%计算实际选择的步长
lamuda =0.999;
for i=1:4
    for j=1:4
        SC(i,j)=R(i,j)/(S(i)*S(j)).^.5;
        HC(i,j)=P(i,j)/(Q(i)*S(j)).^.5;
    end
end
for i=1:4
    for j=1:4
        Dij(i,j)=max([abs(SC(i,j)),abs(HC(i,j)),abs(HC(j,i))]);
    end    
end
Dij=Dij-eye(4);
for i=1:4
Di(i)=max(Dij(i,:));
end
D=max(Di);
if D>0.25
    if k<500
        buchang=0.014;
    else
        buchang=0.004+0.01*D;
    end
elseif D>0.05
    buchang=zeros(4,4);
    for ii=1:4
        if Di(ii)>0.2
            buchang(ii,ii)=0.016;
        elseif Di(ii)>0.05
            buchang(ii,ii)=0.00275+0.05*(Di(ii)-0.05)^0.7;
        else
            buchang(ii,ii)=0.135*Di(ii)^1.3;
        end
    end
else
    for i=1:4
        for j=1:4
            if i~=j
                 buchang(i,j)=0.135*Dij(i,j)^1.3;
            else
                 buchang(i,j)=0.135*Di(i)^1.3;
            end
        end
    end
end



end

