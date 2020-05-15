function subP = updateInformation(subP, Indis, Alg, subSpecies)
%UPDATESUBP �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
if nargin == 4                   % update information from sub=species
    subP.pop(Indis, :) = subSpecies.pop;
    subP.fit(Indis) = subSpecies.fit;
    if strcmp(Alg, 'PSO')
        subP.v(Indis,:) = subSpecies.v;
        subP.cur(Indis,:) = subSpecies.cur;
        subP.curfit(Indis) = subSpecies.curfit;
    end
end

if nargin == 3                  % sort the subP
    subP.pop = subP.pop(Indis, :);
    subP.fit = subP.fit(Indis);
    if strcmp(Alg, 'PSO')
        subP.v = subP.v(Indis, :);
        subP.cur = subP.cur(Indis, :);
        subP.curfit = subP.curfit(Indis);
    end
end

end
