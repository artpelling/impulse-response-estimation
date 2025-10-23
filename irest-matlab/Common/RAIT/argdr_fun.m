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