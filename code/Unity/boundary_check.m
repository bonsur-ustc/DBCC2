function x=boundary_check(x,lbounds,ubounds,dimIndex)
% reflect the individuals over the bounds
lower_bound = lbounds(dimIndex);
upper_bound = ubounds(dimIndex);
length_bound = upper_bound - lower_bound;
x=(x<lower_bound).*(lower_bound+rem((lower_bound-x), length_bound))+(x>=lower_bound).*x;
x=(x>upper_bound).*(upper_bound-rem((x-upper_bound), length_bound))+(x<=upper_bound).*x;

end
