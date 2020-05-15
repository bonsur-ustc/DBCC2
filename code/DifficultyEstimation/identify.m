function relation = identify(group, src_index, des_index, relation, LBounds, LBoundsFitVal, UBounds, func_num)


des_group = [];

for i = 1 : length(des_index)
    des_group = [des_group, group{des_index(i)}];
end


is_interactive = identify_interaction(group{src_index}, des_group, LBounds, LBoundsFitVal, UBounds, func_num);

if is_interactive == 1
    if length(des_index) == 1
        relation(des_index) = src_index;
    else 
        len = length(des_index);
        mid = floor(len/2);
        
        relation = identify(group, src_index, des_index(1:mid), relation, LBounds, LBoundsFitVal, UBounds,func_num);
        
        relation = identify(group, src_index, des_index(mid+1:end), relation, LBounds, LBoundsFitVal, UBounds,func_num);
        
    end
end
end
