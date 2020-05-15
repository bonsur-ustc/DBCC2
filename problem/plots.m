
global initial_flag;
initial_flag = 0;
figure(1);
[X,Y] = meshgrid(-100:100);
Z = zeros(size(X));
i = 1;
for x = -100:100
    j = 1;
    for y = -100:100
        Z(i,j) = func([x y],16,0);
        j = j + 1;
    end
    i = i + 1;
end
 surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');


% ========================================================================
% The  below are changed from CEC'2013 test problems. 
% Please refer to http://www3.ntu.edu.sg/home/EPNSugan/ for original files.
% ==========================================================================
%==============================================================================
%1. Composition Function
%==============================================================================
function fit = CF(x)
global initial_flag
persistent func_num func o sigma lambda bias M


[ps,D] = size(x);
func_num = 5;
lb = -100; ub = 100;
if initial_flag==0
	load Newdata/optima.mat % saved the predefined optima
	if length( o(1,:) ) >= D
		o = o(:,1:D);
	else
		o = lb + (ub - lb) * rand(func_num,D);
	end
	initial_flag=1;
	func.f1 = str2func('FGriewank');
	func.f2 = str2func('FRastrigin');
	func.f3 = str2func('FWeierstrass');
	func.f4 = str2func('FSchwefel');
	func.f5 = str2func('FAckley');
	func.f6 = str2func('FEF8F2');
    func.f7 = str2func('FZakharov');
    func.f8 = str2func('FSphere');
    func.f9 = str2func('FSqure');
    func.f10 = str2func('FSqure');
	bias = zeros(1,func_num);
	sigma = [10,20,30,40,50,60,70,80,90,100];
	lambda = [0.1; 10; 10; 0.1; 100; 2.5; 1e-3;  2.5; 10; 1];
	lambda = repmat(lambda,1,D);
	%for i = 1:func_num
	%	eval(['M.M' int2str(i) '= diag(ones(1,D));']);
	%end
    if D == 2
        load Newdata/CF_M_D2.mat
    end
end
fit = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M);
end


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
it = 0;
res = 0;
for i = 1:func_num
	oo = repmat(o(i,:),ps,1);
	eval(['f = feval(func.f' int2str(i) ',((x-oo)./repmat(lambda(i,:),ps,1))*M.M' int2str(i) ');']);
	eval(['f1 = feval(func.f' int2str(i) ',(x1./lambda(i,:))*M.M' int2str(i) ');']);
	fit1 = 2000 .* f ./ f1;
	res = res + weight(:,i) .* ( fit1 + bias(i) );
end
res = -res;
end


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

%--------------------------------------------------------------------------
% FSchwefel
%--------------------------------------------------------------------------
function f = FSchwefel(x)
[ps,D] = size(x);
f = zeros(ps,1);
for i = 1 : D
    f = f + x(:,i).*sin(sqrt(abs(x(:,i))));
end
f = 418.9829*D*ones(ps,1) - f;
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

%------------------------------------------------------------------------------
% FEF8F2 Function
%------------------------------------------------------------------------------
function f = FEF8F2(x)
[ps,D] = size(x);
f = 0;
for i = 1:(D-1)
	f = f + F8F2( x(:,[i,i+1])+1 );     % (1,...,1) is minimum
end
f = f + F8F2( x(:,[D,1]) +1 );           % (1,...,1) is minimum
end

%------------------------------------------------------------------------------
% F8F2 Function
%------------------------------------------------------------------------------
function f = F8F2(x) 
f2 = 100.*(x(:,1).^2-x(:,2)).^2+(1-x(:,1)).^2; 
f = 1+f2.^2./4000-cos(f2); 
end

%---------------------------------------------------------------------------
% FZakharov Function
%---------------------------------------------------------------------------
function f = FZakharov(x)
[ps,D] = size(x);
sum1 = sum(x.*x,2);
I = 1:D;
I = I';
sum2 = (0.5*x*I).^2;
sum3 = (0.5*x*I).^4;
f = sum1 + sum2 + sum3;
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

% =========================================================================
% generate rotation matrix
% ========================================================================
function M = RotMatrix(D)

A = normrnd(0,1, D,D);
P = orth(A);   % orthogonal matrix
%B = normrnd(0,1,D,D);
%Q = orth(B);   % orthogonal matrix

%u = rand(1, D);
%u = 2.^((u-min(u))./(max(u)-min(u)));
%D = diag(u);


%M = P * D * Q; 
end
