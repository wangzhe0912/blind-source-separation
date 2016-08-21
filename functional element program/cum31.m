function  aa=cum31(a,b)
% a为行向量
a_1 = a.^3;
a_2 = a_1.*b ;
a_3 = a.*b;
a_4 = a.^2;
ea_2 = mean(a_2);
ea_3 = mean(a_4);
ea_4 = mean(a_3);
aa = ea_2 - 3 * ea_3*ea_4;