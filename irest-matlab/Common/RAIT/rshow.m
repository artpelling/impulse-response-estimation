function rshow(s1,s2)

% RSHOW - Visualizes the given function, system or systems.
%
% Usage: 
%     rshow(s)
%     rshow(s1,s2)
%
% Input parameters:
%     s,s1,s2 : complex matrices with rows as elements of a function system
%               (a system with one element is just a plain a function)
%
% Output:
%     Plots of the real and imaginary parts of the elements of the function
%     system.
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

if nargin == 0
    error('Please provide at least on parameter!');
elseif nargin == 1
    [n1,m1] = size(s1);
elseif nargin == 2
    [n1,m1] = size(s1);
    [n2,m2] = size(s2);
    if n1 ~= n2 || m1 ~= m2
        error('The matrices should be of equal size!');
    end
else
    error('Please provide one or two matrices!');
end

x = linspace(0,2*pi,m1+1);
x = x(1:m1);

for i = 1:n1
    switch nargin
        case 1
            subplot(n1,1,i);
            plot(x,real(s1(i,:)),'r',x,imag(s1(i,:)),'b','LineWidth',1);
        case 2
            subplot(n1,2,2*i-1);
            plot(x,real(s1(i,:)),'r',x,imag(s1(i,:)),'b','LineWidth',1);

            subplot(n1,2,2*i);
            plot(x,real(s2(i,:)),'r',x,imag(s2(i,:)),'b','LineWidth',1);
    end
end
