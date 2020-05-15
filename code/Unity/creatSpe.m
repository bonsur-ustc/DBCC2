function subSpecies = creatSpe(subP, Indis, Alg)
%CREATSPE 此处显示有关此函数的摘要
%   此处显示详细说明
% extract the information from subP
subSpecies.pop = subP.pop(Indis,:);
subSpecies.fit = subP.fit(Indis);
if strcmp(Alg, 'PSO')
    subSpecies.v = subP.v(Indis,:);
    subSpecies.cur = subP.cur(Indis, :);
    subSpecies.curfit = subP.curfit(Indis);
end

end

