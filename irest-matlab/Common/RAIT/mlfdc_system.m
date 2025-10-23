function mlf = mlfdc_system(mpoles, eps)

% MLFDC_SYSTEM - Generates the discrete modified basic rational system.
%
% Usage: 
%     mlf = mlfdc_system(mpoles,eps)
%
% Input parameters:
%     mpoles : poles of the discrete modified basic rational system
%     eps    : accuracy of the discretization on the unit disc 
%
% Output parameters:
%     mlf : the elements of the discrete modified basic rational system 
%           at the uniform sampling points as row vectors
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

if max(abs(mpoles)) >= 1
    disp('Bad poles!');
    return;
end
if nargin < 3
    eps = 1e-6;
end

m=length(mpoles);
mlf = zeros(m,m+1); 
t = discretize_dc(mpoles, eps);

[spoles, multi] = multiplicity(mpoles);

for j = 1:length(multi)
    for k=1:multi(j);
        col=sum(multi(1:j-1))+k;
        mlf(col,:)=exp(1i*t).^(k-1)./(1-conj(spoles(j)).*exp(1i*t)).^k;
    end    
end
