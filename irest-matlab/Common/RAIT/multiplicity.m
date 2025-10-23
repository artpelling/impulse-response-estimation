function [spoles, mult] = multiplicity(mpoles)

% MULTIPLICITY - Returns the multiplicity of all elements of 
%                the vector 'mpoles'.
%
% Usage: 
%     [spoles,mult] = multiplicity(mpoles)
%
% Input parameters:
%     mpoles : poles with arbitrary multiplicities
%
% Output parameters:
%     spoles : vector of the poles that contains a pole only once
%     mult   : mult(i) refers to the multiplicity of the ith pole
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

poles =  mpoles;
mult = zeros(size(mpoles));
spoles = zeros(size(mpoles));
top=0;

while ~isempty(poles)
    ind = ones(1,length(poles));
    top = top + 1;
    for i = 1:length(poles)
        if poles(1) == poles(i)
            ind(i) = i;
            mult(top) = mult(top) + 1;
        end
    end
    spoles(top)=poles(1);
    poles(ind)=[];
end
mult = mult(1:top);
spoles = spoles(1:top);