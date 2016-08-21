function [ SIR1,SIR2 ] = functionSIR( s,y )
%计算SIR信干比

s1=s(1,:);s2=s(2,:);y1=y(1,:);y2=y(2,:);
s11=2*(s1-min(s1))/(max(s1)-min(s1))-1;
s22=2*(s2-min(s2))/(max(s2)-min(s2))-1;
is1=2*(max(y1)-y1)/(max(y1)-min(y1))-1;
is2=2*(max(y2)-y2)/(max(y2)-min(y2))-1;
S=sum(is1.*s11,2);
tt=sum(s11.*s11,2);
ta=S*s11/tt;
Si=sum(is2.*s22,2);
tti=sum(s22.*s22,2);
tai=Si*s22/tti;
Rss=[sum(s11.*s11,2),sum(s11.*s22,2);sum(s22.*s11,2),sum(s22.*s22,2)];
gy=inv(Rss);
jy=[sum(is1.*s11,2),sum(is1.*s22,2)]';
ty=(gy*jy)';
ss=[s11;s22];
re=ty*ss-ta;
fr=[sum(is2.*s11,2),sum(is2.*s22,2)]';
rew=(gy*fr)';
tq=rew*ss-tai;
SIR1=10*log10((sum(ta.*ta,2))/(sum(re.*re,2)))
SIR2=10*log10((sum(tai.*tai,2))/(sum(tq.*tq,2)))


end

