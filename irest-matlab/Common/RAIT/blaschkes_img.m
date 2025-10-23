function [B_abs,B_arg,B]=blaschkes_img(path, a, show)

% BLASCHKES_IMG - Transforms an image by applying the Blaschke function 
%                 defined by poles a.
%
% Usage: 
%     b = blaschkes_img(path,a,show)
%
% Input parameters:
%     path : path of the input image 
%     a    : parameters of the Blaschke function
%     show : displays the transformed images
%
% Output parameters:
%     B_abs : absolute values of the Blaschke function 
%     B_arg : arguments of the Blaschke function 
%     B     : transformed image
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

M=imread(path);
M=M(:,:,1);
B=0.5.*ones(size(M));
B_arg=0.5.*ones(size(M));
B_abs=0.5.*ones(size(M));
[m,n]=size(M);
blaschke=@(z) (z-a)./(1-conj(a).*z);
%blaschke_inv=@(z) (z+a)./(1-conj(-a).*z);
u=size(M,1);
half=round(u/2);
bound=half-10;

for i=1:1:m
    for j=1:1:n
        if sqrt((j-half)^2+(i-half)^2)<bound
            z=(j-half)/bound + (i-half)*1i/bound;
            b=blaschke(z);
            I=round(bound*real(b)+half);
            J=round(bound*imag(b)+half);
            B(i,j)=M(J,I);
            B_abs(i,j)=abs(blaschke(z));
            t=atan(imag(blaschke(z))/real(blaschke(z)));
            B_arg(i,j)=t;            
        else
            B(i,j)=NaN;
%            B_arg(i,j)=NaN;
%            B_abs(i,j)=NaN;            
        end
    end
end

if show
    %Image%
    figure;    
    imshow(uint8(B));
    hold on;
    half=round(length(B)/2);
    plot(real(a)*(half-10)+half,imag(a)*(half-10)+half,'r.','MarkerSize',30);
%    title('Transformed image');

    %Absolute value of Blaschke function%
    figure;        
    imshow(B_abs);
    colormap(jet);    
%    title('Absolute value of Blaschke function');       

    %Argument function%
    figure;   
    imshow(B_arg);
    colormap(jet);
%    title('Argument function');
end