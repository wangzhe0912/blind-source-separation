function hh = RIR2( RR,SS )
%计算混合系数
%   此处显示详细说明
c = 340;                    % Sound velocity (m/s)
fs = 8000;                 % Sample frequency (samples/s)
L = [7 5 2.75];                % Room dimensions [x y z] (m)
beta = 0.2;                 % Reverberation time (s)
n = fs*beta;                   % Number of samples
N=size(RR,2);
M=size(SS,2);
height=1.5;

receivers=[3.96,1,height;
           4.04,1,height;
           3.88,1,height;
           4.12,1,height;];
sources=[3.8500,1.2598,height;
         3.7181,1.1026,height;
         4.1026,1.2819,height;
         4.2598,1.1500,height;];

     
for i=1:N 
    switch RR(i)
        case 1
            R(i,:)=receivers(1,:);
        case 2
            R(i,:)=receivers(2,:);
        case 3
            R(i,:)=receivers(3,:);
        case 4
            R(i,:)=receivers(4,:);
    end
end

for i=1:M 
    switch SS(i)
        case 1
            S(i,:)=sources(1,:);
        case 2
            S(i,:)=sources(2,:);
        case 3
            S(i,:)=sources(3,:);
        case 4
            S(i,:)=sources(4,:);
    end
end

for i=1:N
    for j=1:M
        r=R(i,:);
        s=S(j,:);
        hh(i,j,:) = rir_generator(c, fs, r, s, L, beta, n);
    end
end

end