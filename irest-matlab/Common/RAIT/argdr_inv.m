function t = argdr_inv(a,b,epsi)

% ARG_INV - Inverse images by the argument function of a Blaschke product.
%
% Usage: 
%     t = arg_inv(a,b)
%     t = arg_inv(a,b,epsi)
%     
%
% Input parameters:
%     a       : parameters of the Blaschke product
%     b       : values in [-pi,pi) whose inverse image is needed
%     epsi : required precision for the inverses (optional, default 1e-4)
%
% Output parameters:
%     t : inverse images by the argument function of the points in t
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

if nargin == 2
    epsi = 1e-4;
end

na = size(a,1);
nb = size(b,1);
if na ~= 1 || nb ~= 1
    error('Parameters should be row vectors!');
end
if max(abs(a)) >= 1
    error('Bad poles!');
end


if length(a) == 1
    t = arg_inv_one(a,b);
else
    t = arg_inv_all(a,b,epsi);
end

% -------------------------------------------------------------------------

function t = arg_inv_one(a,b)

r = abs(a);
fi = angle(a);
mu = (1+r)/(1-r);

gamma = 2 * atan((1/mu)*tan(fi/2));

t = 2*atan((1/mu)*tan((b-gamma)/2)) + fi;
t = mod(t+pi,2*pi)-pi; % move it in [-pi,pi)


% -------------------------------------------------------------------------

function x = arg_inv_all(a,b,epsi)

% uses the bisection method with an enhanced order of calculation

n = length(b);
s = bisection_order(n) + 1;
x = zeros(1,n+1);
for i = 1:n+1

    if i == 1
        v1 = -pi;
        v2 = pi;
        fv1 = -pi; % fv1 <= y
        fv2 = pi;  % fv2 >= y
    elseif i == 2
        x(n+1) = x(1) + 2*pi; %x(s(2,1))
        continue;
    else
        v1 = x(s(i,2));
        v2 = x(s(i,3));
        fv1 = (argdr_fun(a,v1)-v1/2)./length(a);
        fv2 = (argdr_fun(a,v2)-v2/2)./length(a);
    end
    
    ba = b(s(i,1));
    if fv1 == ba
        x(s(i,1)) = v1;
        continue;
    elseif fv2 == ba
        x(s(i,1)) = v2;
        continue;
    else
        xa = (v1+v2)/2;
        fvk = (argdr_fun(a,xa)-xa/2)./length(a);
        while abs(fvk-ba) > epsi
            if fvk == ba
                x(s(i,1)) = xa;
                return
            elseif fvk < ba
                v1 = xa;
            else
                v2 = xa;
            end
            xa = (v1+v2)/2;
            fvk = (argdr_fun(a,xa)-xa/2)./length(a);
        end
        x(s(i,1)) = xa;
    end
end
x = x(1:n);

% -------------------------------------------------------------------------

%Same function as 'arg_fun', but it is continuous on IR. 
function b = argdr_fun(a,t)

b = zeros(1,length(t));
b=zeros(1,length(t));
for j=1:1:length(t)
    for i = 1:1:length(a);
        bs=arg_fun(a(i),t(j));
        b(j)=b(j)+bs+2*pi*floor((t(j)+pi)/(2*pi));
    end
end
% -------------------------------------------------------------------------