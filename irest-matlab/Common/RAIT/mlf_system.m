function mlf = mlf_system(len, mpoles)

% Generates the modified lf system defined by 'poles'.

% MLF_SYSTEM - Generates the modified basic rational function system 
%              defined by mpoles.
%
% Usage: 
%     mlf = mlf_system_(len,mpoles)
%
% Input parameters:
%     len    : number of points in case of uniform sampling 
%     mpoles : poles of the modified basic rational system
%
% Output parameters:
%     mlf : the elements of the modified basic rational function system 
%           at the uniform sampling points as row vectors
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

[np,mp] = size(mpoles);
if np ~= 1 || len < 2
    disp('Wrong parameters!');
    return;
end
if max(abs(mpoles)) >= 1
    disp('Bad poles!');
    return;
end

mlf = zeros(mp,len);
t = linspace(-pi,pi,len+1);
t = t(1:len);
z = exp(1i*t);

[spoles, multi] = multiplicity(mpoles);

for j = 1:length(multi)
    for k=1:multi(j);
        col=sum(multi(1:j-1))+k;
        mlf(col,:)=z.^(k-1)./(1-conj(spoles(j)).*z).^k;
    end    
end