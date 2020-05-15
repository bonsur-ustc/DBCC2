function subSpecies = creatSpe(subP, Indis, Alg)
%CREATSPE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% extract the information from subP
subSpecies.pop = subP.pop(Indis,:);
subSpecies.fit = subP.fit(Indis);
if strcmp(Alg, 'PSO')
    subSpecies.v = subP.v(Indis,:);
    subSpecies.cur = subP.cur(Indis, :);
    subSpecies.curfit = subP.curfit(Indis);
end

end

