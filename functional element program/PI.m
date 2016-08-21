function [ E ] = PI( P )


guimo=size(P);
n=guimo(1);
sumi1=0;
for i=1:n
    sumj1=0;
    for j=1:n
        sumj1=sumj1+abs(P(i,j))/max(abs(P(i,:)));
    end
    k=sumj1-1;
    sumi1=sumi1+k;
end
part1=sumi1/n;

sumj2=0;
for j=1:n
    sumi2=0;
    for i=1:n
        sumi2=sumi2+abs(P(i,j))/max(abs(P(:,j)));
    end
    k=sumi2-1;
    sumj2=sumj2+k;
end
part2=sumj2/n;


E=part1+part2;
end

