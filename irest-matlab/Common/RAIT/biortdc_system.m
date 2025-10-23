function bts = biortdc_system(mpoles, eps)

% BIORTDC_SYSTEM - Generates the discrete biorthogonal system.
%
% Usage: 
%     bts = biortdc_system(mpoles,eps)
%
% Input parameters:
%     mpoles : poles of the discrete biorthogonal system
%     eps    : accuracy of the discretization on the unit disc 
%
% Output parameters:
%     bts : the elements of the discrete biorthogonal system at the uniform 
%           sampling points as row vectors
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

[np,mp] = size(mpoles);
if max(abs(mpoles)) >= 1
    error('Bad poles!');
    return;
end
if nargin < 3
    eps = 1e-6;
end

m=length(mpoles);
bts = zeros(m,m+1); 
t = discretize_dc(mpoles, eps);

[spoles, multi] = multiplicity(mpoles);

for j = 1:length(multi)
    for k=1:multi(j);
        col=sum(multi(1:j-1))+k;
        bts(col,:)=pszi(j,k,spoles,multi,exp(1i*t));
    end    
end

% -------------------------------------------------------------------------
% The following functions and notations are used from :
%
% S. Fridli, F. Schipp, “Biorthogonal systems to rational functions,” 
% Annales Univ. Sci. Budapest., Sect. Comp, vol. 35, no. 1, 
% pp. 95–105, 2011.


function v=pszi(l,j,poles,multi,z)

% Compute the values of the biorthognal polynomial at z related to the lth
% pole with j multiplicity.

n=length(poles);
v=zeros(1,length(z));
Do=Domega(multi(l)-j,l,poles,multi,poles(l));

for s=0:1:multi(l)-j
    v=v+Do(s+1)./factorial(s).*(z-poles(l)).^s;
end

v=v.*Omega(l,poles,multi,z)./Omega(l,poles,multi,poles(l)).*(z-poles(l)).^(j-1);

% -------------------------------------------------------------------------

function v=Omega(l,poles,multi,z)

% Computes the values of the Omega base functions related to the
% biorthogonal system.

n=length(poles);
v=ones(1,length(z));
v=v ./ ( 1 - conj(poles(l)).*z ).^multi(l);
% Blaschke-function
B=@(z,a) (z-a) ./ (1-conj(a).*z); 

for i=1:1:l-1
    v=v.*B(z,poles(i)).^multi(i);
end
for i=l+1:1:n
    v=v.*B(z,poles(i)).^multi(i);
end

% -------------------------------------------------------------------------

function Do=Domega(s,l,poles,multi,z)

% Computes sth derivative of the omega function. 
% The first row of array 'Do' contains omega.%
% The rth derivative is stored in Do(r+1,:). 

n=length(poles);
Do=zeros(s+1,length(z));
Do(1,:)=Omega(l,poles,multi,poles(l))./Omega(l,poles,multi,z);

for i=1:1:s
    for j=1:1:i
        Do(i+1,:)=Do(i+1,:)+nchoosek(i-1,j-1).*Do(j,:).*ro(i-j,l,poles,multi,z);
    end
end

% -------------------------------------------------------------------------

function v=ro(s,l,poles,multi,z)

% Computes sth derivative of the auxiliary function to Domega.

n=length(multi);
v=ones(1,length(z));
v=multi(l) ./ (z-(1/conj(poles(l)))).^(s+1);

for i=1:1:l-1
    v=v-multi(i).*( (1./(z-poles(i))).^(s+1) - (1./(z-(1/conj(poles(i))))).^(s+1) );
end
for i=l+1:1:n
    v=v-multi(i).*( (1./(z-poles(i))).^(s+1) - (1./(z-(1/conj(poles(i))))).^(s+1) );
end

v=v.*(-1).^s.*factorial(s);

