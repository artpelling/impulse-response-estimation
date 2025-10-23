%Construct the test signal.%
t_1=-pi/4   ;
t_2=pi/4 + pi/6;
t_3=-pi/2 - pi/5;
mpoles=[0.85*exp(1i*t_1), 0.85*exp(1i*t_1),  0.73*exp(1i*(t_2)), 0.73*exp(1i*(t_3))];
c1=(2.1+0.15*1i);
c2=(2.1+0.3*1i);
c3=(1.1-0.3*1i);
co=[0 c1 c2 c3];
lfs = lf_system(1000, mpoles);
sig_real=real(co*lfs);
sig_imag=addimag(sig_real);
sig=sig_imag;

lfs = lf_system(1000, mpoles);
disp('% lf_system -------- OK');

sig_lf = lf_generate(1000, mpoles, co);
disp('% lf_generate ------ OK');
%plot(real(sig_lf),'b');

mlfs = mlf_system(1000, mpoles);
disp('% mlf_system ------- OK');

[mlf_co, mlf_err] = mlf_coeffs(sig, mpoles);
disp('% mlf_coeffs ------- OK');

sig_mlf = mlf_generate(1000, mpoles, mlf_co);
disp('% mlf_generate ----- OK');
% plot(real(sig_mlf),'r');

mlfdc = mlfdc_system(mpoles);
disp('% mlfdc_system ----- OK');

[mlfdc_co, mlfdc_err] = mlfdc_coeffs(sig, mpoles);
disp('% mlfdc_coeffs ----- OK');

sig_mlfdc = mlfdc_generate(1000, mpoles, mlfdc_co);
disp('% mlfdc_generate --- OK');
%plot(real(sig_mlfdc),'r');

mt = mt_system(1000, mpoles);
disp('% mt_system -------- OK');

[mt_co, mt_err] = mt_coeffs(sig, mpoles);
disp('% mt_coeffs -------- OK');

sig_mt = mt_generate(1000, mpoles, mt_co);
disp('% mt_generate ------ OK');
%plot(1:1:1000, real(sig), 'g', 1:1:1000, real(sig_mt), 'r');

mtdc = mtdc_system(mpoles);
disp('% mtdc_system ------ OK');

[mtdc_co, mtdc_err] = mtdc_coeffs(sig, mpoles);
disp('% mtdc_coeffs ------ OK');

sig_mtdc = mtdc_generate(1000, mpoles, mtdc_co);
disp('% mtdc_generate ---- OK');
%plot(1:1:1000, real(sig), 'g', 1:1:1000, real(sig_mtdc), 'r');

mtdr = mtdr_system(mpoles);
disp('% mtdr_system ------ OK');

[mtdr_co_re, mtdr_co_im, mtdr_err] = mtdr_coeffs(real(sig), mpoles);
disp('% mtdr_coeffs ------ OK');

sig_mtdr = mtdr_generate(1000, mpoles, mtdr_co_re, mtdr_co_im);
disp('% mtdr_generate ---- OK');
%plot(1:1:1000, real(sig), 'g', 1:1:1000, sig_mtdr, 'r');

bts = biort_system(1000, mpoles);
disp('% biort_system ----- OK');

[bts_co, bts_err] = biort_coeffs(sig, mpoles);
disp('% biort_coeffs ----- OK');

sig_bts = biort_generate(1000, mpoles, bts_co);
disp('% biort_generate --- OK');
%plot(real(sig_bts),'g')

btsdc = biortdc_system(mpoles);
disp('% biortdc_system --- OK');

[btsdc_co, btsdc_err] = biortdc_coeffs(sig, mpoles);
disp('% biortdc_coeffs --- OK');

sig_btsdc = biort_generate(1000, mpoles, btsdc_co);
disp('% biortdc_generate - OK');
%plot(real(sig_btsdc),'g')

mt_co=coeff_conv(1000, mpoles, co, 'lf', 'mt');
disp('% coeff_conv ------- OK');
%sig_mt = mt_generate(1000, mpoles, mt_co);
%plot(1:1:1000, real(sig_lf), 'g', 1:1:1000, real(sig_mt), 'r--');

biortdc_co=coeffd_conv(mpoles, mtdc_co, 'mtdc', 'biortdc');
disp('% coeffd_conv ------ OK');
%sig_biortdc = biortdc_generate(1000, mpoles, biortdc_co);
%plot(1:1:1000, real(sig_lf), 'g', 1:1:1000, real(sig_biortdc), 'r--');

mult=[1 2 1];  
period=1;
[p,c] = simplex_mt(sig_real, mult, period, [0.3 0.3 0.4 0.1 0.5 0.2], 0, 1e-6);
disp('% simplex_mt ------- OK');
% %mpoles=periodize_poles(multiply_poles(p,mult),period)
% %plot(1:1:1000, real(sig_lf), 'g', 1:1:1000, real(mt_generate(1000,mpoles,c)), 'r--');

mult=[1 2 1];  
period=1;
[p,c] = simplex_mtdc(sig_real, mult, period, [0.3 0.3 0.4 0.0 0.5 0.2], 0, 1e-3);
disp('% simplex_mtdc ----- OK');
% %mpoles=periodize_poles(multiply_poles(p,mult),period)
% %plot(1:1:1000, real(sig_lf), 'g', 1:1:1000, real(mtdc_generate(1000,mpoles,c)), 'r--');

mult=[1 2 1];  
period=2;
[p,cu,cv] = simplex_mtdr(sig_real, mult, period, [-0.5 -0.4 0.0 -0.0 0.3 0.4], 0, 1e-3);
disp('% simplex_mtdr ----- OK');
% mpoles=periodize_poles(multiply_poles(p,mult),period)
% plot(1:1:1000, real(sig_lf), 'g', 1:1:1000, real(mtdr_generate(1000,mpoles,cu,cv)), 'r--');

mult=[1 2 1];  
period=1;
[p,c] = simplex_biort(sig_real, mult, period, [1 0 1 0 1 0], 0, 1e-10);
disp('% simplex_biort ---- OK');
% mpoles=periodize_poles(multiply_poles(p,mult),period)
% plot(1:1:1000, real(sig_lf), 'g', 1:1:1000, real(biort_generate(1000,mpoles,c)), 'r--');

mult=[1 4 1];  
period=1;
[p,c] = simplex_biortdc(sig_real, mult, period, [0.6 -0.6 0 0 0.2 0.6], 0, 1e-6);
disp('% simplex_biortdc -- OK');
%mpoles=periodize_poles(multiply_poles(p,mult),period)
%plot(1:1:1000, real(sig_lf), 'g', 1:1:1000, real(biortdc_generate(1000,mpoles,c)), 'r--');