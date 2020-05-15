function fit = func(x, fun_num, iseval,dimIndex, base)
%*********************************************************************
% last modified: 2019/8/19
% separable many-modal problems
%*********************************************************************
% iseval decides whether evaluations are countable
global evals;
global initial_flag;

[pSize, D] = size(x);
if nargin == 5 && length(dimIndex) < size(base,2) 
    % base = lbound
    tmp_x = repmat(base, pSize, 1);
    tmp_x(:,dimIndex) = x;
    x = tmp_x;
    D = size(base, 2);
end
P = 1 : D; % dimension


switch(fun_num)
    case 1
        % nonseparable
        subPro = [2 2 2 2 2];
        func_handles = {@CF1,@CF3,@CF5,@CF7,@CF9};
        eval(['load data/IM_M1.mat']);
        x = x*M1;
    case 2
        % nonseparable
        subPro = [2 2 2 2 2];
        func_handles = {@CF2,@CF4,@CF6,@CF8,@CF10};
        eval(['load data/IM_M2.mat']);
        x = x*M2;  
    case 3
        subPro = [2,2,2,2,2];
        func_handles = {@CF3,@CF4,@CF6,@CF6,@CF7};
    case 4
        subPro = [2,2,2,2,2];
        func_handles = {@CF3,@CF4,@CF5,@CF6,@CF8};
    case 5
        subPro = [3, 3, 2, 2];
        func_handles = {@CF1,@CF2,@CF3,@CF5};
    case 6
        subPro = [2,3,2,3];
        func_handles = {@CF5,@CF6,@CF6,@CF8};
    case 7 
        subPro = [2,3,2,3];
        func_handles = {@CF6,@CF7,@CF8,@CF9};
    case 8
        subPro = [3,2,2,3];
        func_handles ={@CF1,@CF5,@CF6,@CF7};
    case 9
        subPro = [2,3,2,3];
        func_handles = {@CF2,@CF3,@CF5,@CF8}; 
    case 10
        subPro = [2,3,2,3];
        func_handles = {@CF2,@CF5,@CF5,@CF9}; 
    case 11
        subPro = [3,2,2,3];
        func_handles ={@CF1,@CF4,@CF8,@CF9};
    case 12 
        subPro = [2,3,2,3];
        func_handles = {@CF1,@CF2,@CF8,@CF10};
    case 13 
        subPro = [2,3,2,3];
        func_handles = {@CF2,@CF8,@CF9,@CF10};
    case 14
        subPro = [4, 2, 2, 2];
        func_handles = {@CF1, @CF2, @CF3, @CF8};
    case 15
        subPro = [4,2,2,2];
        func_handles = {@CF2,@CF8,@CF9,@CF10};
end
%%
fit = 0;
ldim = 0;
persistent init_count;
for i = 1 : length(subPro)
    if initial_flag == 0
        init_count = 1;                %load the information 
    end
    dims = P(ldim+1:ldim+subPro(i));
    ldim = ldim + subPro(i);
    fit = fit + func_handles{i}(x(:,dims));
    if init_count < length(subPro)
        initial_flag = 0;
        init_count = init_count + 1;
    end
end
MAX = 1;
fit = MAX * fit;
if iseval
    evals = evals + pSize;  %count the fitness evaluation
end

end

%%
% ========================================================================
% The  below are changed from CEC'2013 test problems. 
% Please refer to http://www3.ntu.edu.sg/home/EPNSugan/ for original files.
% ==========================================================================
%==============================================================================
%1. Composition Function 1, n = 5
%==============================================================================
function fit = CF1(x)
global initial_flag
persistent func_num func o sigma lambda bias M


[ps,D] = size(x);
func_num = 5;
lb = -100; ub = 100;
if initial_flag==0
	eval(['load data/optima.mat']); % saved the predefined optima
	if length( o(1,:) ) >= D
		o = o(:,1:D);
	else
		o = lb + (ub - lb) * rand(func_num,D);
	end
	initial_flag=1;
	func.f1 = str2func('FGriewank');
	func.f2 = str2func('FRastrigin');
	func.f3 = str2func('FWeierstrass');
	func.f4 = str2func('FSphere');
	func.f5 = str2func('FSqure');
	bias = zeros(1,func_num); 
	sigma = [10,20,60,20,10];
	lambda = [1; 10; 20; 1; 1];
	lambda = repmat(lambda,1,D);
    c = ones(1,func_num);
	%for i = 1:func_num
	%	eval(['M.M' int2str(i) '= diag(ones(1,D));']);
	%end
    if D == 2
        eval(['load data/CF1_M_D2.mat']);
    elseif D == 3
        eval(['load data/CF1_M_D3.mat']);
    elseif D == 4
        eval(['load data/CF1_M_D4.mat']);
    elseif D == 5
        eval(['load data/CF1_M_D5.mat']);
    elseif D == 10 
        eval(['load data/CF1_M_D10.mat']);
    elseif D == 20
        eval(['load data/CF1_M_D20.mat']);
    else
        for i = 1:func_num
            eval(['M.M' int2str(i) '=RotMatrixCondition(D,c(i));']);
        end
    end
end
fit = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M);
end

function fit = CF2(x)
global initial_flag;
persistent func_num func o sigma lambda bias M
[ps,D] = size(x);
func_num = 5;
lb = -100; ub = 100;

if initial_flag == 0
    eval(['load data/optima.mat']); % saved the predefined optima
    if length( o(1,:) ) >= D
        o = o(:,1:D);
    else
        o = lb + (ub - lb) * rand(func_num,D);
    end
    initial_flag = 1;
    func.f1 = str2func('FAckley');
    func.f2 = str2func('FPeriodic');
    func.f3 = str2func('FSphere');
    func.f4 = str2func('FSalomon');
    func.f5 = str2func('FWeierstrass');
    bias = zeros(1,func_num);
    sigma = [30,20,10,10,40];
    lambda = [10;10;10;10;20];
    lambda = repmat(lambda,1,D);
    c = ones(1,func_num);
    if D == 2
        eval(['load data/CF2_M_D2.mat']);
    elseif D == 3
        eval(['load data/CF2_M_D3.mat']);
    elseif D == 4
        eval(['load data/CF1_M_D4.mat']);
    elseif D == 5
        eval(['load data/CF2_M_D5.mat']);
    elseif D == 10
        eval(['load data/CF2_M_D10.mat']);
    elseif D == 20
        eval(['load data/CF2_M_D20.mat']);
    else
        for i = 1:func_num
            eval(['M.M' int2str(i) '= RotMatrixCondition(D,c(i));']);
        end
    end
end
fit = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M);
end


function fit = CF3(x)
global initial_flag;
persistent func_num func o sigma lambda bias M
[ps,D] = size(x);
func_num = 2;
lb = -100; ub = 100;

if initial_flag == 0
    eval(['load data/optima.mat']); % saved the predefined optima
    if length( o(1,:) ) >= D
        o = o(:,1:D);
    else
        o = lb + (ub - lb) * rand(func_num,D);
    end
    initial_flag = 1;
    func.f1 = str2func('FAckley');
    func.f2 = str2func('FAckley');
    bias = zeros(1,func_num);
    sigma = [10,20];
    lambda = [200;200];
    lambda = repmat(lambda,1,D);
    c = ones(1,func_num);
    if D == 2
        eval(['load data/CF3_M_D2.mat']);
%         for i = 1:func_num
%             eval(['M.M' int2str(i) '= diag(ones(1,D));']);
%         end
    elseif D == 3
        eval(['load data/CF3_M_D3.mat']);
    elseif D == 5
        eval(['load data/CF3_M_D5.mat']);
    elseif D == 10
        eval(['load data/CF3_M_D10.mat']);
    elseif D == 20
        eval(['load data/CF3_M_D20.mat']);
    else
        for i = 1:func_num
            eval(['M.M' int2str(i) '= RotMatrixCondition(D,c(i));']);
        end
    end
end
fit = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M);
end

function fit = CF4(x)
global initial_flag;
persistent func_num func o sigma lambda bias M
[ps,D] = size(x);
func_num = 2;
lb = -100; ub = 100;

if initial_flag == 0
    eval(['load data/optima.mat']); % saved the predefined optima
    if length( o(1,:) ) >= D
        o = o(:,1:D);
    else
        o = lb + (ub - lb) * rand(func_num,D);
    end
    initial_flag = 1;
    func.f1 = str2func('FSphere');
    func.f2 = str2func('FAckley');
    bias = zeros(1,func_num);
    sigma = [100,100];
    lambda = [10;40];
    lambda = repmat(lambda,1,D);
    c = ones(1,func_num);
    if D == 2
        eval(['load data/CF4_M_D2.mat']);
%         for i = 1:func_num
%             eval(['M.M' int2str(i) '= diag(ones(1,D));']);
%         end
    elseif D == 3
        eval(['load data/CF4_M_D3.mat']);
    elseif D == 5
        eval(['load data/CF4_M_D5.mat']);
    elseif D == 10
        eval(['load data/CF4_M_D10.mat']);
    elseif D == 20
        eval(['load data/CF4_M_D20.mat']);
    else
        for i = 1:func_num
            eval(['M.M' int2str(i) '= RotMatrixCondition(D,c(i));']);
        end
    end
end
fit = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M);
end


function fit = CF5(x)
global initial_flag;
persistent func_num func o sigma lambda bias M
[ps,D] = size(x);
func_num = 3;
lb = -100; ub = 100;

if initial_flag == 0
    eval(['load data/optima.mat']); % saved the predefined optima
    if length( o(1,:) ) >= D
        o = o(:,1:D);
    else
        o = lb + (ub - lb) * rand(func_num,D);
    end
    initial_flag = 1;
    func.f1 = str2func('FSphere');
    func.f2 = str2func('FSphere');
    func.f3 = str2func('FGriewank');


    bias = zeros(1,func_num);
    sigma = [15,15,15];
    lambda = [1;1;10];
    lambda = repmat(lambda,1,D);
    c = ones(1,func_num);
    if D == 2
        eval(['load data/CF5_M_D2.mat']);
    elseif D == 3
        eval(['load data/CF5_M_D3.mat']);
    elseif D == 5
        eval(['load data/CF5_M_D5.mat']);
    elseif D == 10
        eval(['load data/CF5_M_D10.mat']);
    elseif D == 20
        eval(['load data/CF5_M_D20.mat']);
    else
        for i = 1:func_num
            eval(['M.M' int2str(i) '= RotMatrixCondition(D,c(i));']);
        end
    end
end
fit = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M);
end


function fit = CF6(x)
global initial_flag;
persistent func_num func o sigma lambda bias M
[ps,D] = size(x);
func_num = 3;
lb = -100; ub = 100;

if initial_flag == 0
    eval(['load data/optima.mat']); % saved the predefined optima
    if length( o(1,:) ) >= D
        o = o(:,1:D);
    else
        o = lb + (ub - lb) * rand(func_num,D);
    end
    initial_flag = 1;
    func.f1 = str2func('FSphere');
    func.f2 = str2func('FPeriodic');
    func.f3 = str2func('FWeierstrass');
    bias = zeros(1,func_num);
    sigma = [30,20,40];
    lambda = [1;10;60];
    lambda = repmat(lambda,1,D);
    c = ones(1,func_num);
    if D == 2
        eval(['load data/CF6_M_D2.mat']);
    elseif D == 3
        eval(['load data/CF6_M_D3.mat']);
    elseif D == 5
        eval(['load data/CF6_M_D5.mat']);
    elseif D == 10
        eval(['load data/CF6_M_D10.mat']);
    elseif D == 20
        eval(['load data/CF6_M_D20.mat']);
    else
        for i = 1:func_num
            eval(['M.M' int2str(i) '= RotMatrixCondition(D,c(i));']);
        end
    end
end
fit = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M);
end

function fit = CF7(x)
global initial_flag;
persistent func_num func o sigma lambda bias M
[ps,D] = size(x);
func_num = 3;
lb = -100; ub = 100;

if initial_flag == 0
    eval(['load data/optima.mat']); % saved the predefined optima
    if length( o(1,:) ) >= D
        o = o(:,1:D);
    else
        o = lb + (ub - lb) * rand(func_num,D);
    end
    initial_flag = 1;
    func.f1 = str2func('FGriewank');
    func.f2 = str2func('FRastrigin');
    func.f3 = str2func('FPeriodic');
    bias = zeros(1,func_num);
    sigma = [60,80,80];
    lambda = [1;10;10];
    lambda = repmat(lambda,1,D);
    c = ones(1,func_num);
    if D == 2
        eval(['load data/CF7_M_D2.mat']);
    elseif D == 3
        eval(['load data/CF7_M_D3.mat']);
    elseif D == 5
        eval(['load data/CF7_M_D5.mat']);
    elseif D == 10
        eval(['load data/CF7_M_D10.mat']);
    elseif D == 20
        eval(['load data/CF7_M_D20.mat']);
    else
        for i = 1:func_num
            eval(['M.M' int2str(i) '= RotMatrixCondition(D,c(i));']);
        end
    end
end
fit = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M);
end

function fit = CF8(x)
global initial_flag;
persistent func_num func o sigma lambda bias M
[ps,D] = size(x);
func_num = 4;
lb = -100; ub = 100;

if initial_flag == 0
    eval(['load data/optima.mat']); % saved the predefined optima
    if length( o(1,:) ) >= D
        o = o(:,1:D);
    else
        o = lb + (ub - lb) * rand(func_num,D);
    end
    initial_flag = 1;
    func.f1 = str2func('FGriewank');
    func.f2 = str2func('FRastrigin');
    func.f3 = str2func('FSphere');
    func.f4 = str2func('FSqure');
    bias = zeros(1,func_num);
    sigma = [60,20,20,50];
    lambda = [10;20;10;30];
    lambda = repmat(lambda,1,D);
    c = ones(1,func_num);
    if D == 2
        eval(['load data/CF8_M_D2.mat']);
    elseif D == 3
        eval(['load data/CF8_M_D3.mat']);
    elseif D == 5
        eval(['load data/CF8_M_D5.mat']);
    elseif D == 10
        eval(['load data/CF8_M_D10.mat']);
    elseif D == 20
        eval(['load data/CF8_M_D20.mat']);
    else
        for i = 1:func_num
            eval(['M.M' int2str(i) '= RotMatrixCondition(D,c(i));']);
        end
    end
end
fit = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M);
end


function fit = CF9(x)
global initial_flag;
persistent func_num func o sigma lambda bias M
[ps,D] = size(x);
func_num = 4;
lb = -100; ub = 100;

if initial_flag == 0
    eval(['load data/optima.mat']); % saved the predefined optima
    if length( o(1,:) ) >= D
        o = o(:,1:D);
    else
        o = lb + (ub - lb) * rand(func_num,D);
    end
    initial_flag = 1;
    func.f1 = str2func('FGriewank');
    func.f2 = str2func('FGriewank');
    func.f3 = str2func('FWeierstrass');
    func.f4 = str2func('FAckley');
    bias = zeros(1,func_num);
    sigma = [50,60,70,80];
    lambda = [1/4;1/4;80;30];
    lambda = repmat(lambda,1,D);
    c = ones(1,func_num);
    if D == 2
        eval(['load data/CF8_M_D2.mat']);
    elseif D == 3
        eval(['load data/CF8_M_D3.mat']);
    elseif D == 5
        eval(['load data/CF8_M_D5.mat']);
    elseif D == 10
        eval(['load data/CF8_M_D10.mat']);
    elseif D == 20
        eval(['load data/CF8_M_D20.mat']);
    else
        for i = 1:func_num
            eval(['M.M' int2str(i) '= RotMatrixCondition(D,c(i));']);
        end
    end
end
fit = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M);
end



function fit = CF10(x)
global initial_flag;
persistent func_num func o sigma lambda bias M
[ps,D] = size(x);
func_num = 4;
lb = -100; ub = 100;

if initial_flag == 0
    eval(['load data/optima.mat']); % saved the predefined optima
    if length( o(1,:) ) >= D
        o = o(:,1:D);
    else
        o = lb + (ub - lb) * rand(func_num,D);
    end
    initial_flag = 1;
    func.f1 = str2func('FAckley');
    func.f2 = str2func('FPeriodic');
    func.f3 = str2func('FWeierstrass');
    func.f4 = str2func('FWeierstrass');
    bias = zeros(1,func_num);
    sigma = [100,20,40,80];
    lambda = [10;10;50;30];
    lambda = repmat(lambda,1,D);
    c = ones(1,func_num);
    if D == 2
        eval(['load data/CF8_M_D2.mat']);
    elseif D == 3
        eval(['load data/CF8_M_D3.mat']);
    elseif D == 5
        eval(['load data/CF8_M_D5.mat']);
    elseif D == 10
        eval(['load data/CF8_M_D10.mat']);
    elseif D == 20
        eval(['load data/CF8_M_D20.mat']);
    else
        for i = 1:func_num
            eval(['M.M' int2str(i) '= RotMatrixCondition(D,c(i));']);
        end
    end
end
fit = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M);
end


%%
%==============================================================================
% Hybrid Composition General Framework
%==============================================================================
function res = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M)
[ps,D] = size(x);
x1 = 100*ones(1,D);
for i = 1:func_num
	oo = repmat( o(i,:), ps, 1 );
	weight(:,i) = exp( -sum( (x-oo).^2, 2 )./2./( D * sigma(i)^2 ) );
end

[tmp,tmpid] = sort(weight,2);
for i = 1:ps
	weight(i,:) = (weight(i,:)==tmp(i,func_num)) .* weight(i,:) + (weight(i,:)~=tmp(i,func_num)) .* (weight(i,:).*(1-tmp(i,func_num).^10));
end
if sum(weight,2) == 0
	weight = weight + 1;
end
weight = weight ./ repmat( sum(weight,2), 1, func_num );
res = 0;
for i = 1:func_num
	oo = repmat(o(i,:),ps,1);
	eval(['f = feval(func.f' int2str(i) ',((x-oo)./repmat(lambda(i,:),ps,1))*M.M' int2str(i) ');']);
	eval(['f1 = feval(func.f' int2str(i) ',(x1./lambda(i,:))*M.M' int2str(i) ');']);
	fit1 = 1000 .*f ./ abs(f1);
	res = res + weight(:,i) .* ( fit1 + bias(i) );
end
res = -res;
end

%%
%==============================================================================
% Basic Functions
%==============================================================================

%------------------------------------------------------------------------------
% Sphere Function
%------------------------------------------------------------------------------
function f = FSphere(x)
%Please notice there is no use to rotate a sphere function, with rotation
%here just for a similar structure as other functions and easy programming
[ps,D] = size(x);
f = sum( x.^2, 2);
end

%------------------------------------------------------------------------------
% Griewank's Function
%------------------------------------------------------------------------------
function f = FGriewank(x)
[ps,D] = size(x);
f = 1;
for i = 1:D
	f = f.*cos( x(:,i)./sqrt(i) );
end
f = sum( x.^2, 2)./4000 - f + 1;
end

%------------------------------------------------------------------------------
% Rastrigin's Function
%------------------------------------------------------------------------------
function f = FRastrigin(x)
[ps,D] = size(x);
f = sum( x.^2-10.*cos( 2.*pi.*x )+10, 2 );
end

%------------------------------------------------------------------------------
% Weierstrass Function
%------------------------------------------------------------------------------
function f = FWeierstrass(x)
[ps,D] = size(x);
x = x+0.5;
a = 0.5;
b = 3;
kmax = 20;
c1(1:kmax+1) = a.^(0:kmax);
c2(1:kmax+1) = 2*pi*b.^(0:kmax);
f = 0;
c = -w(0.5,c1,c2);
for i = 1:D
	f = f + w( x(:,i)', c1, c2 );
end
f = f + c*D;
end

function y = w(x,c1,c2)
y = zeros(length(x),1);
for k = 1:length(x)
	y(k) = sum( c1.*cos( c2.*x(:,k) ) );
end
end

%----------------------------------------------------------------------------
% FAckley
%----------------------------------------------------------------------------
function f = FAckley(x)
[ps,D] = size(x);
a = 20;
b = 0.2;
c = 2*pi;
sum1 = sum(x.*x,2)./D;
sum2 = sum(cos(c*x),2)./D;
f = -a*exp(-b*sqrt(sum1))-exp(sum2) + a + exp(1);
end

%--------------------------------------------------------------------------
% Sum Squares function
%-------------------------------------------------------------------------
function f = FSqure(x)
[ps,D] = size(x);
f = zeros(ps,1);
tmp = x.*x;
I = 1 :D;
I = I'; 
f = tmp*I;
end

function f = FPeriodic(x) % modified minimum = 0
sum1 = sin(x).^2;
sum2 = sum(x.^2,2);
f = sum(sum1,2)-0.1*exp(-sum2) + 0.1;  %
end

function f = FSalomon(x)
x2 = x .^ 2;
sumx2 = sum(x2, 2);
sqrtsx2 = sqrt(sumx2);

f = 1 - cos(2 .* pi .* sqrtsx2) + (0.1 * sqrtsx2);
end


% ====================================================================
% Please refer to http://www3.ntu.edu.sg/home/EPNSugan/ for original files.
% ===================================================================
function [q,r] = LocalGramSchmidt (A)
% computes the QR factorization of $A$ via
% classical Gram Schmid 

[n,m] = size(A); 
q = A;    
for j=1:m
    for i=1:j-1 
        r(i,j) = q(:,j)'*q(:,i);
    end
    for i=1:j-1   
      q(:,j) = q(:,j) -  r(i,j)*q(:,i);
    end
    t =  norm(q(:,j),2 ) ;
    q(:,j) = q(:,j) / t ;
    r(j,j) = t  ;
end 
end

%------------------------------------------------------------------------------
% Generates a D-dimensional rotation matrix with predifined Condition Number (c)
% Please refer to http://www3.ntu.edu.sg/home/EPNSugan/ for original files.
%------------------------------------------------------------------------------
function M = RotMatrixCondition(D,c)

% A random normal matrix
A = normrnd(0,1,D,D);

% P Orthogonal matrix
P = LocalGramSchmidt(A);

% A random normal matrix
A = normrnd(0,1,D,D);

% Q Orthogonal matrix
Q = LocalGramSchmidt(A);

% Make a Diagonal matrix D with prespecified Condition Number
u = rand(1,D);
D = c .^ ( (u-min(u))./(max(u)-min(u)) );
D = diag(D);

% M rotation matrix with Condition Number c
M = P * D * Q;
end

% =========================================================================
% generate rotation matrix
% ========================================================================
function M = RotMatrix(D)

A = normrnd(0,1, D,D);
M = orth(A);   % orthogonal matrix 
end
