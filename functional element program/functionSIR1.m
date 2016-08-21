function [ SIR ] = functionSIR1( s,y )
%计算SIR信干比
     [N,k]=size(s);
     Rss=zeros(N,N);
%      s(1,:)=2*(s(1,:)-min(s(1,:)))/(max(s(1,:))-min(s(1,:)))-1;
%      s(2,:)=2*(s(2,:)-min(s(2,:)))/(max(s(2,:))-min(s(2,:)))-1;
%      y(1,:)=2*(max(y(1,:))-y(1,:))/(max(y(1,:))-min(y(1,:)))-1;
%      y(2,:)=2*(max(y(2,:))-y(2,:))/(max(y(2,:))-min(y(2,:)))-1;
     for i=1:N
         for j=1:N
            Rss(i,j)=dot(s(i,:),s(j,:));
         end
     end
     Rss_1=inv(Rss);
     s_target=zeros(N,k);e_interf=zeros(N,k);
     for jj=1:N
         s_target(jj,:)=dot(y(jj,:),s(jj,:))*s(jj,:)/norm(s(jj,:))^2;
         for ii=1:N
             row_vector(ii)=dot(y(jj,:),s(ii,:));
         end
         sum(jj,:)=(Rss_1*row_vector')'*s;
         e_interf(jj,:)=sum(jj,:)-s_target(jj,:);
         SIR(jj)=10*log10(norm(s_target(jj,:))^2/norm(e_interf(jj,:))^2);
     end
     SIR
end

