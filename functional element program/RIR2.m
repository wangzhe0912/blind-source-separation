function hh = RIR2( RR,SS )
%计算混合系数
%   此处显示详细说明
c = 340;                    % Sound velocity (m/s)
fs = 8000;                 % Sample frequency (samples/s)
L = [7 5 2.75];                % Room dimensions [x y z] (m)
beta = 0.2;                 % Reverberation time (s)
n = 4096;                   % Number of samples
N=size(RR,2);
M=size(SS,2);
height=1.5;

receivers=[3.91,1,height;
           3.97,1,height;
           4.03,1,height;
           4.09,1,height];
sources=[2.52,1.26,height;
         3.75,2.30,height;
         4,2.5,height;
         4,3,height;
         4.25,2.3,height;
         4.96,2.15,height;
         5.28,2.53,height;
         5.41,1.51,height;
         5.88,1.68,height;
         5.48,1.26,height];

     
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
        case 5
            S(i,:)=sources(5,:);
        case 6
            S(i,:)=sources(6,:);
        case 7
            S(i,:)=sources(7,:);
        case 8
            S(i,:)=sources(8,:);
        case 9
            S(i,:)=sources(9,:);
        case 10
            S(i,:)=sources(10,:);
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