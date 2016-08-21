function yy=leijiliang(y1,y2,L)
b=0;
%for j=1:L
  %  if j==1
      % t = buyiwei(j);
   % else
   %     d=j-1;
     %   t=yiwei(d);
   % end
  %      bb=filter(t,1,y1);
   %     ym1=bb;
for k = 1:L
    if k==1
       a = buyiwei(k);
    else
        c=k-1;
        a=yiwei(c);
    end
        dd=filter(a,1,y2);
        ym2=dd;
    y_k(1,:)=y1;
    y_k(2,:)=ym2;
    %y_k = white(y_k);                    % °×»¯¡ª¡ªºÜ¹Ø¼ü
    yw1 = y_k(1,:);
    yw2 = y_k(2,:);
    b=b+abs(cum31(yw1,yw2))+abs(cum22(yw1,yw2));
    %b=b+abs(cum31(ym1,ym2)+cum22(ym1,ym2));
   %end
end
       yy=1/b;
        % yy=b;
            