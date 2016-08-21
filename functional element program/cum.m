function  aa=cum(a)
% a为行向量
yy1=a(1,:);
yy2=a(2,:);
yy1_2=yy1.^2;
yy2_2=yy2.^2;
e1=mean(yy1_2.*yy2_2);
e2=mean(yy1_2)*mean(yy2_2);
e3=2*((mean(yy1.*yy2))^2);
aa=1/(e1-e2-e3);



%a_2 = a.^2;
%a_4 = a.^4;
%ea_2 = mean(a_2);
%ea_4 = mean(a_4);
%aa = ea_4 - 3 * ea_2 ^2;