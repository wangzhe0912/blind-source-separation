function  aa=cum22(a,b)
% a为行向量
a_2 = a.^2;
a_4 = b.^2;
a_1=a_2.*a_4;
a_3=a.*b;
ea_1 = mean(a_1);
ea_2 = mean(a_2);
ea_3 = mean(a_4);
ea_4 = mean(a_3);
aa = ea_1 - ea_2*ea_3 - 2 * (ea_4^2);