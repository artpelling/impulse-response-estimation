t_1=-pi/4;
t_2=pi/4 + pi/6;
t_3=-pi/2 - pi/5;
mpoles=[0.85*exp(1i*t_1) 0.85*exp(1i*t_1),  0.73*exp(1i*(t_2)), 0.73*exp(1i*(t_3))];
c1=(2.1+0.15*1i);
c2=(2.1+0.3*1i);
c3=(1.1-0.3*1i);
co=[0 c1 c2 c3];
lf_co=co;
lfs = lf_system(1000, mpoles);
sig_real=real(co*lfs);
sig_imag=addimag(sig_real);
sig=sig_imag;

mt_co=coeff_conv(1000, mpoles, co, 'lf', 'mt');
sig_mt = mt_generate(1000, mpoles, mt_co);
sig_lf = lf_generate(1000, mpoles, co);
subplot(3,2,1);
hold on;
plot(real(sig_lf),'g','LineWidth',2);
plot(real(sig_mt),'r');
[~,err1] = mt_coeffs(sig_lf,mpoles);
legend(num2str(err1));
title('lf2mt')

mlf_co=coeff_conv(1000, mpoles, lf_co, 'lf', 'mlf');
sig_lf = lf_generate(1000, mpoles, lf_co);
sig_mlf = mlf_generate(1000, mpoles, mlf_co);
subplot(3,2,2);
hold on;
plot(real(sig_lf),'g','LineWidth',2);
plot(real(sig_mlf),'r');
[~,err1] = mlf_coeffs(sig_lf,mpoles);
legend(num2str(err1));
title('lf2mlf')

biort_co=coeff_conv(1000, mpoles, lf_co, 'lf', 'biort');
sig_lf = lf_generate(1000, mpoles, lf_co);
sig_biort = biort_generate(1000, mpoles, biort_co);
subplot(3,2,3);
hold on;
plot(real(sig_lf),'g','LineWidth',2);
plot(real(sig_biort),'r');
[~,err1] = biort_coeffs(sig_lf,mpoles);
legend(num2str(err1));
title('lf2biort')

biort_co=coeff_conv(1000, mpoles, mlf_co, 'mlf', 'biort');
sig_mlf = mlf_generate(1000, mpoles, mlf_co);
sig_biort = biort_generate(1000, mpoles, biort_co);
subplot(3,2,4);
hold on;
plot(real(sig_mlf),'g','LineWidth',2);
plot(real(sig_biort),'r');
[~,err1] = biort_coeffs(sig_mlf,mpoles);
legend(num2str(err1));
title('mlf2biort')

biort_co=coeff_conv(1000, mpoles, mt_co, 'mt', 'biort');
sig_mt = mt_generate(1000, mpoles, mt_co);
sig_biort = biort_generate(1000, mpoles, biort_co);
subplot(3,2,5);
hold on;
plot(real(sig_mt),'g','LineWidth',2);
plot(real(sig_biort),'r');
[~,err1] = biort_coeffs(sig_mt,mpoles);
legend(num2str(err1));
title('mt2biort');

mlf_co=coeff_conv(1000, mpoles, mt_co, 'mt', 'mlf');
sig_mt = mt_generate(1000, mpoles, mt_co);
sig_mlf = mlf_generate(1000, mpoles, mlf_co);
subplot(3,2,6);
hold on;
plot(real(sig_mt),'g','LineWidth',2);
plot(real(sig_mlf),'r');
[~,err1] = mlf_coeffs(sig_mt,mpoles);
legend(num2str(err1));
title('mt2mlf')

% -------------------------------------------------------------------------
figure;
mtdc_co=mtdc_coeffs(real(sig_mt), mpoles, 1e-6);
biortdc_co=coeffd_conv(mpoles, mtdc_co, 'mtdc', 'biortdc',1e-6);
sig_mt = mt_generate(1000, mpoles, mtdc_co);
sig_biort = biort_generate(1000, mpoles, biortdc_co);
subplot(3,1,1);
hold on;
plot(real(sig_mt),'g','LineWidth',2);
plot(real(sig_biort),'r');
[~,err1] = mtdc_coeffs(sig_mt,mpoles);
legend(num2str(err1));
title('mtdc2biortdc')

mtdc_co=mtdc_coeffs(real(sig_mt), mpoles, 1e-6);
mlfdc_co=coeffd_conv(mpoles, mtdc_co, 'mtdc', 'mlfdc',1e-6);
sig_mt = mt_generate(1000, mpoles, mtdc_co);
sig_mlf = mlf_generate(1000, mpoles, mlfdc_co);
subplot(3,1,2);
hold on;
plot(real(sig_mt),'g','LineWidth',2);
plot(real(sig_mlf),'r');
[~,err1] = mlfdc_coeffs(sig_mt,mpoles);
legend(num2str(err1));
title('mtdc2mlfdc')

[mlfdc_co, err]=mlfdc_coeffs(real(sig_mt), mpoles, 1e-6);
biortdc_co=coeffd_conv(mpoles, mlfdc_co, 'mlfdc', 'biortdc', 1e-6);
sig_mlf = mlf_generate(1000, mpoles, mlfdc_co);
sig_biort = biort_generate(1000, mpoles, biortdc_co);
subplot(3,1,3);
hold on;
plot(real(sig_mlf),'g','LineWidth',2);
plot(real(sig_biort),'r');
[c,err1] = biortdc_coeffs(sig_mlf,mpoles);
legend(num2str(err1));
title('mlfdc2biortdc');

