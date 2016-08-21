function [ SIR ] = SIR( y,s,L )
%根据y和s计算分离算法的性能
[n,k]=size(y);
for i=1:n
    for j=1:n
        [s_target,e_interf,e_artif]=bss_decomp_filt(y(j,:),i,s,L);
        [SDR(i,j),SIR(i,j),SAR(i,j)]=bss_crit(s_target,e_interf,e_artif);
    end
end
SIR

end

