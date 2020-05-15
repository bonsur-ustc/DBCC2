function finalgroups = fast_separating(groups, lb, lbFitVal, ub,func_num)

group_num = size(groups, 2);
finalgroups={};
while group_num>0
    if group_num ==1
        finalgroups=[finalgroups,groups{1}];
        groups(1)=[];
        group_num=group_num-1;
        continue;
    end
    
    % first variable
    src_index=1;
    
    % rest variables
    des_index=2:group_num;
    
    % relation(i) = j if  group i  is interactive with group j
    relation=zeros(1,group_num);
    relation = identify(groups, src_index, des_index, relation, lb, lbFitVal, ub,func_num);
    
    index=find(relation==1);
    len=length(index);
    
    % add interactive variavles to first group
    for i=1:len
        groups{1}=[groups{1},groups{index(i)}];
    end
    
    % delete interactive group
    groups(index)=[];
    group_num=group_num-length(index);
    
    % inseparable variables are record in final groups
    if len==0
        finalgroups=[finalgroups,groups{1}];
        groups(1)=[];
        group_num=group_num-1;
    end
end
end
