function vi = addimag(v)

% ADDIMAG - Calculates the imaginary part of v using FFT.
%
% Usage:
%     vi = addimag(v)
%
% Input paraameters:
%     v : a vector with real elements
%
% Output parameters:
%     vi : a complex vector with appropriate imaginary part
%          (to be in Hardy space)
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??


if size(v,1) ~= 1
    error('Wrong vector!');
end
if max(imag(v)) ~= 0
    error('The vector is not real!');
end

vf = fft(v);
vif = mt_arrange(vf);
vi = ifft(vif);


% -------------------------------------------------------------------------

function ta = mt_arrange(t)

% Rearrage FFT(v) so that lots of zeros appear on the right side of the FFT.

mt = size(t,2);
ta = zeros(size(t));
ta(1) = t(1);
for i = 2:mt/2
    ta(i) = t(i) + conj(t(mt+2-i));
    ta(mt+2-i) = 0;
end

