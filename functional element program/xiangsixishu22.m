function ans1 = xiangsixishu22(s,y)
bb=[s(1,:)' y(1,:)'];
dd=[s(2,:)' y(2,:)'];
bb1=[s(1,:)' y(2,:)'];
dd1=[s(2,:)' y(1,:)'];
%qq=[s(3,:)' y(3,:)'];
format long
cc=corrcoef(bb) ; % 源信号与恢复的信号的相关系数
ee=corrcoef(dd);
cc2=corrcoef(bb1) ; % 源信号与恢复的信号的相关系数
ee2=corrcoef(dd1);
ans1=[cc(1,2) cc2(1,2); ee2(1,2)  ee(1,2)];