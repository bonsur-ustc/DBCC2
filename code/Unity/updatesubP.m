function subP = updatesubP(subP, subSpecies, index, Alg)
%UPDATESUBP �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
subP.pop(index, :) = subSpecies.pop;
subP.fit(index) = subSpecies.fit;
if strcmp(Alg, 'PSO')
    subP.v(index,:) = subSpecies.v;
    subP.cur(index,:) = subSpecies.cur;
    subP.curfit(index) = subSpecies.curfit;
end
end

