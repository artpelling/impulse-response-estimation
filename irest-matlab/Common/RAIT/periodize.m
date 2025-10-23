function s=periodize(v, alpha, draw)

% PERIODIZE - Calculates the periodized extension of a signal.
%
% Usage: 
%     s = periodize(v,alpha,draw)
%
% Input parameters:
%     v     : an arbitrary vector  
%     alpha : alpha parameter of the tukey window (tukeywin)
%     draw  : optional logical value to display the periodized signal
%
% Output parameters:
%     s  : periodized signal 
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

if size(v,1)~=1
    error('Signal must be row vector!')
end
if nargin<3
    draw=false;
end

[v_max,v_ind]=max(v);
gap=(max(v(1),v(end))+min(v(1),v(end)))/2;
v=v-gap;
n=length(v);
N=ceil(n*alpha);
%approximating derivatives%
smooth=sgolayfilt(v,2,11);
fx=remNEP(smooth);
[~,dy]=curvatures(1:1:length(v),smooth);
left=sum(dy(1:fx(2)));
right=sum(dy(fx(end-1):end));
s=[ones(1,N).*v(1),v,ones(1,N).*v(end)];

if sign(left)==sign(right)
    endslope=(v(end)-v(1))/n; 
    if sign(left)==sign(endslope)
        %taking into account the significance of
        %the first derivatives at the end points%
        trs=abs(max(v(1),v(end))-min(v(1),v(end)))/2;
        if abs(left)<abs(right)
            s=s-sign(right)*(abs(v(end))+trs);
        elseif abs(left)>abs(right)
            s=s+sign(left)*(abs(v(1))+trs);        
        end
    else
        avg=v(1)+v(end)/2;           
        s=s-avg;
    end    
    
elseif sign(left)>0
    %positive first derivatives%
    if v(1)<0 || v(end)<0
        s=s-min(v(1),v(end));
    end
else 
    %negative first derivatives%
    if sign(v(1))>0 || sign(v(end))>0
        s=s-max(v(1),v(end));
    end
end
tukey=tukeywin(length(s),(2*N+2)/length(s));
s=s.*tukey';
s=s+(v_max-s(N+v_ind));

if draw
    plot(s,'r');
    hold on;
    plot(N+1:1:N+length(v),s(N+1:1:N+length(v)),'g');    
    hold off;    
end

% -------------------------------------------------------------------------

function [x,y] = remNEP(data)

% Remove Non Extreme Points from the given data array.

n=length(data);
j=1;
x=zeros(n,1);
y=zeros(n,1);
x(1)=1;
y(j)=data(j);
for i=2:1:n-1
    %if data(i) min or max extreme point
    if (y(j)< data(i) && data(i)>data(i+1)) || (y(j)> data(i) && data(i)<data(i+1))
        y(j+1)=data(i);
        x(j+1)=i;
        j=j+1;
    end    
end
y(j+1)=data(n);
x(j+1)=n;
x=x(1:j+1);
y=y(1:j+1);
  
function [k,dy,ddy]=curvatures(x,y)

% Approximating discrete curvatures.

n=length(y);
dy=zeros(n-2,1);
ddy=zeros(n-2,1);
%Setting the curvatures to zero at the endpoints.%
k=zeros(n-2,1);
for i=2:1:n-1
    %Bézier-approximation%
    t=x(i);
    ddy(i-1)=y(i-1).*(-2/(x(i+1)-x(i-1))).*(-1/(x(i+1)-x(i-1))) - 2.*(4.*y(i)-y(i-1)-y(i+1))/(x(i+1)-x(i-1)).^2 + y(i+1).*(2/(x(i+1)-x(i-1)).^2);
    dy(i-1)=2.*y(i-1).*((x(i+1)-t)/(x(i+1)-x(i-1))).*(-1/(x(i+1)-x(i-1))) + (4.*y(i)-y(i-1)-y(i+1)).*(x(i-1)+x(i+1)-2.*x(i))/(x(i+1)-x(i-1)).^2 + 2.*y(i+1).*((t-x(i-1))/(x(i+1)-x(i-1))).*(1/(x(i+1)-x(i-1)));
    k(i-1)=ddy(i-1)./((1+dy(i-1).^2).^(3/2));
end
