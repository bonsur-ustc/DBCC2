% ******************************************************************************
% * Version: 1.0
% * Last modified on: 21 January, 2013 
% * Developers: Michael G. Epitropakis, Xiaodong Li.
% *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
% *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
% * ****************************************************************************
function [count, finalseeds] = count_goptima(pop, nfunc, accuracy)

% pop: NP, D
[NP, D] = size(pop);

% information of test functions
rho = 0.5*ones(1,15);
fgoptima = zeros(1,15);
nopt = [360,480,108,144,150,108,144,135,120,180,160,400,320,200,320];


% evaluate pop
% fpop = zeros(1,NP);
% for i=1:NP
% 	fpop(i) = mixed(pop(i,:), nfunc);
% end
% fpoptmp = fpop;
% fpop = fit; % modified
fpop = func(pop,nfunc,false);
fpoptmp = fpop;

% descent sorting
[B, IX] = sort(fpoptmp,'descend');

% Sort population based on its fitness values
% do not change the current populatio population. Work on cpop/cpopfits
cpop = pop(IX,:);
cpopfits = fpop(IX);

%get seeds
seeds = [];
seedsidx = [];

for i=1:NP
	found=0;
	[sNP,sD] = size(seeds);
	for j=1:sNP
		% Calculate distance from seeds
		dist = sqrt( sum( (seeds(j,:)-cpop(i,:)).^2,2) );
		% If the Euclidean distance is less than the radius
		if (dist <= rho(nfunc))
			found = 1;
			break;
		end
	end
	% If it is not similar to any other seed, then it is a new seed
	if (found == 0)
		seeds = [seeds;cpop(i,:)];
		seedsidx = [seedsidx; i];
	end
end

% Based on the accuracy: check which seeds are global optimizers
count = 0; finalseeds = [];
seedsfit = cpopfits(seedsidx);
[ idx ] = find(abs(seedsfit - fgoptima(nfunc))<=accuracy);
if (length(idx) > nopt(nfunc) )
	idx = idx(1:nopt(nfunc));
end
count = length(idx);
finalseeds = seeds(idx,:);
