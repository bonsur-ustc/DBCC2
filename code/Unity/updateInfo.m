function subP = updateInfo(subP, Index, Alg)
%UPDATEP �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
subP.pop = subP.pop(Index, :);
subP.fit = subP.fit(Index);
if strcmp(Alg, 'PSO')
    subP.v = subP.v(Index, :);
    subP.cur = subP.cur(Index, :);
    subP.curfit = subP.curfit(Index);
end
end

