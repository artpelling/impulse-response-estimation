function arg_inv_anim(a,n)

% ARG_INV_ANIM - Shows an animation related to the inverse of an
% equidistant discretization by an argument function
%
% Usage: 
%     arg_inv_anim(a,n)
%
% Input parameters:
%     a : the parameter of a Blaschke function
%     n : number of discretization points
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??


[na,ma] = size(a);
[nn,mn] = size(n);
if na ~= 1 || ma ~= 1 || nn ~= 1 || mn ~= 1
    error('Please provide two numbers!');
end
if abs(a) >= 1
    error('The parameter should be inside the unit disc!');
end

t = linspace(-pi,pi,n+1); % discretization
t = t(1:n);

anim = 32;
part = 2*pi/n/anim;
curr = 0;

for i = 1:anim*n
    
    b = arg_inv(a,t+curr);
    plot(exp(1i*b),'ko');
    hold on;
    plot(real(a),imag(a),'ro');
    hold off;
    axis equal;
    axis([-1.2 1.2 -1.2 1.2]);
    
    curr = curr + part;
    if curr > 2*pi/n
        curr = curr - 2*pi/n;
    end
    pause(0.01);
end
