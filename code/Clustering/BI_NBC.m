function [species,cutfit] = BI_NBC(pop, popfit)
%  NBC divides population into multi sub-popultions
%  For more information please refer to the
%  M. Preuss. "Niching the CMA-ES via nearest-better clustering." In
%  Proceedings of the 12th annual conference companion on Genetic and evolutionary
%  computation (GECCO ¡¯10). ACM, New York, NY, USA, pp. 1711-1718, 2010.

cutfit = mean(popfit);          %the value of cutoff  
cutfit = cutfit - 1e-10; 

normal = find(popfit >= cutfit);
%normal = 1:size(pop,1);                                 % size normal = NP means BI_NBC is the same as NBC 
matdis = pdist2(pop(normal,:),pop(normal,:)); % distance matrix

factor= 2.0;        % fai
n=length(normal);
nbc=zeros(n,3);
nbc(1:n,1)=1:n;
nbc(1,2)=-1;          % nbc(i,2) records the nearest better index of individual i
nbc(1,3)=0;             % nbc(i,3) for distance

for i=2:n
    [u,v]=min(matdis(i,1:i-1)); % u for the minimum distance£¬v for index
    nbc(i,2)=v;
    nbc(i,3)=u;
end

meandis=factor*mean(nbc(2:n,3)); % the mean distance
nbc(nbc(:,3)>meandis,2)=-1;      % if larger than meandis, then cut it
nbc(nbc(:,3)>meandis,3)=0;
seeds=nbc(nbc(:,2)==-1,1);

m=zeros(n,2);
m(1:n,1)=1:n;
for i=1:n                % find the root
    j=nbc(i,2);
    k=j;
    while j~=-1
        k=j;
        j=nbc(j,2);
    end
    if k==-1
        m(i,2)=i;
    else
        m(i,2)=k;
    end
end

species = struct();   % construct the result
for i=1:length(seeds)
    species(i).seed = seeds(i);
    species(i).idx = m(m(:, 2) == seeds(i), 1);
    species(i).len = length(species(i).idx);
end


total = 1 :size(pop,1);
outliers = setdiff(total,normal);
% to do
% outliers = length(normal)+1 : total;

if ~isempty(outliers)   
    dis = pdist2(pop(outliers,:),pop(seeds,:));
    for i = 1 : length(outliers)
        [~,near] = min(dis(i,:));
        species(near).idx = [species(near).idx;outliers(i)];
        species(near).len = species(near).len + 1;
    end
end






