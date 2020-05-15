function is_interactive = identify_interaction(src_group, des_group, LBounds, LBoundsFitVal, UBounds,func_num)
is_interactive = 0;

dim = size(LBounds, 2);
% m = 0.5*(lb + ub);

MBounds = 0.5*(LBounds+UBounds);
 
% ind2(i) = m(i)
ind2=LBounds;
ind2(src_group)=MBounds(src_group);
val2=func(ind2,func_num,true);

% ind3(j) = m(j)
ind3=LBounds;
ind3(des_group)=MBounds(des_group);
val3=func(ind3,func_num,true);

% ind4(i,j) = m(i,j)
ind4=LBounds;
ind4([src_group,des_group])=MBounds([src_group,des_group]);
val4=func(ind4,func_num,true);

% delta(i,j) = abs(delta1 - delta2)
delta1=LBoundsFitVal-val2;
delta2=val3-val4;
muM = eps / 2;
gamma = @(n)((n.*muM)./(1-n.*muM));

% error sup
errub = gamma(dim^0.5) * max([abs(LBoundsFitVal),abs(val2),abs(val3),abs(val4)]);
alpha=abs(delta1-delta2);
if(alpha>errub)
    is_interactive = 1;
end
end


