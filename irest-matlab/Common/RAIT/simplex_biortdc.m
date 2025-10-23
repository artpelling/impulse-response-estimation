function [p,c] = simplex_biortdc(f, mult, period, inits, show, epsilon)

% SIMPLEX_BIORTDC - Gives the poles of the discrete biorthogonal system 
%                   that best fit the approximation of the function 'f'.
%
% Usage: 
%     [p,c] = simplex_biortdc(f,mult,period,inits,show,epsilon)
%
% Input parameters:
%     f       : an arbitrary vector  
%     mult    : multiplicities of the desired poles
%     period  : periodicity of the desired poles
%     inits   : values of the initial simplex
%     show    : optional logical value to display the approximation process 
%     epsilon : accuracy of the process
%
% Output parameters:
%     p : predicted poles of the discrete biorthogonal system with 
%         'mult' multiplicities and 'period' periodicity
%     c : the Fourier coefficients of 'f' with respect to the discrete 
%         biorthogonal system defined by the predicted poles 'p' 
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??


% initial simplex
polenum = length(mult); % pólusok száma
n = 2*polenum; % dimension
ne = n + 1; % number of vertices of the simplex

% simplex:
% vertices1coord1 vertices1coord2 ... vertices1coordn function value1
% vertices2coord1 vertices2coord2 ... vertices2coordn function value2
% ...
% verticesne                      ...                 function valuene

simplex = zeros(ne,ne);
h = .1;
simplex(1,:) = [inits 0];
he          = h * ones(1,1+polenum*2);
for j = 2:ne
    simplex(j,:) = simplex(1,:);
    simplex(j,j-1) = simplex(j,j-1) + he(j-1);
end

alfa = 1; % reflection parameter
beta = 0.5; % narrowing parameter
gamma = 2; % extending parameter
if nargin<6
    epsilon=1e-6;
end

% for display

styles = ['bo'; 'bx'; 'b.'; 'bd';'bp'];


% -------------------------------------------------------------------------


% approximated function

f = addimag(f);

% approximated signal

sample = size(f,2);

t = linspace(-pi,pi,sample+1);
t = t(1:sample);
unit_disc = exp(1i*t);


% values of the initial simplex

for j = 1:ne
    [~,simplex(j,ne)] = biortdc_coeffs(f,periodize_poles(multiply_poles(coords2params(simplex(j,1:n)),mult),period)); % függvényérték az eltérés
end

% iteration

step = 0;
while std(simplex(:,ne)) > epsilon
    step = step + 1;
    
	% ordering
	simplex = sortrows(simplex,ne);
	
	% centroid
	xs = sum(simplex(1:n,1:n)) / n;
	
	% reflected point
	xr = (1 + alfa) * xs - alfa * simplex(n+1,1:n);
	[~,fvxr] = biortdc_coeffs(f,periodize_poles(multiply_poles(coords2params(xr),mult),period)); % fvxr
	
	if fvxr < simplex(1,ne) % best function value   
        % examining the extended point
        xe = gamma * xr + (1 - gamma) * xs;
        [~,fvxe] = biortdc_coeffs(f,periodize_poles(multiply_poles(coords2params(xe),mult),period)); % fvxe
        
        if fvxe < fvxr
            simplex(ne,:) = [xe,fvxe];
        else
            simplex(ne,:) = [xr,fvxr];
        end
	end
	
	if fvxr >= simplex(1,ne) && fvxr <= simplex(n,ne) % 2. közepes függvényérték
        simplex(ne,:) = [xr,fvxr];
	end
	
	if fvxr > simplex(n,ne)
        % examining the narrowed point
        if fvxr >= simplex(ne,ne)
            xc = beta * simplex(ne,1:n) + (1 - beta) * xs;
        else
            xc = beta * xr + (1 - beta) * xs;
        end
        [~,fvxc] = biortdc_coeffs(f,periodize_poles(multiply_poles(coords2params(xc),mult),period)); % fvxc
        
        if fvxc < fvxr && fvxc < simplex(ne,ne)
            simplex(ne,:) = [xc,fvxc];
        else
            for k = 2:ne
                simplex(k,1:n) = (simplex(k,1:n) + simplex(1,1:n))/2;
                [~,simplex(k,ne)] = biortdc_coeffs(f,periodize_poles(multiply_poles(coords2params(simplex(k,1:n)),mult),period));
            end
        end
	end
	
    % display
    
    if show
    
        subplot(1,2,1);
        sz = coords2params_all(simplex(:,1:n));
        c = biortdc_coeffs(f,periodize_poles(multiply_poles(coords2params(simplex(1,1:n)),mult),period));
        fs = biortdc_generate(sample,periodize_poles(multiply_poles(coords2params(simplex(1,1:n)),mult),period),c);
        m=sum(mult)*period;
        z=linspace(-pi,pi,m+1);
        argt=arg_inv(periodize_poles(multiply_poles(coords2params(simplex(1,1:n)),mult),period),z, epsilon);
        smp=subsample(f, argt);
        plot(unit_disc,'k');
        title(['step: ' int2str(step)]);
        hold on;
        for k = 1:polenum
            plot(sz(:,k),styles(k,:));
        end
        % plot(c,'rd');
        % plot(poles,'rd');
		hold off;
        
        subplot(1,2,2);
		plot(t,zeros(1,sample),'k');
		hold on;
        plot(t,real(f),'r');
		% plot(t,real(f),'m',t,imag(f),'c');
		plot(t,real(fs),'b');
        stem(argt,smp,'b');  
        % plot(t,real(fs),'r',t,imag(fs),'b');
		% plot(t,real(fs-f),'g',t,imag(fs-f),'g:');
		hold off;
	    drawnow;
        
        if step == 1
            pause
        end
        if rem(step,10) == 0
            pause(0.01);
        end
        pause(0.01);
    end
    
end


% Output:

simplex = sortrows(simplex,ne);
p = coords2params(simplex(1,1:n));
c = biortdc_coeffs(f,periodize_poles(multiply_poles(p,mult),period));


