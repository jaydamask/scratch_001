function h = make_h_bessel_ct(M1, order, T, dt)

%--------------------------------------------------------------------------------------------------
%
% fcxn:    make_h_bessel_ct.m
% descrip: constructs a continuous-time bessel-filter impulse response 
%          parameterized by /M1/ and sampled on the interval [0: Nwindow-1].
% domain:  continuous 
% inputs:  /M1/       location parameter, the first moment
%          /order/    order of pole
%          /T/        final time
%          /dt/       time interval
% output:  /h/        impulse response on [0:dt:T-dt] of bessel
% author:  JN Damask
%
%--------------------------------------------------------------------------------------------------

% sanity
if (order < 1) || (order > 5),
    disp('order out of bounds:   1 <= order <= 5');
    return
end

% make taxis
taxis = [0: dt: T-dt];

% dictionary of Bessel polynomials
bessel.poly{1} = [1 1];
bessel.poly{2} = [1 3 3];
bessel.poly{3} = [1 6 15 15];
bessel.poly{4} = [1 10 45 105 105];
bessel.poly{5} = [1 15 105 420 945 945];

% compute pfe
bessel.last_coeff = max(bessel.poly{order});
[A, s, k] = residue(bessel.last_coeff, bessel.poly{order});

% add up the partial impulse responses
h = real(A.' * exp((s / M1) * taxis)) / M1;
h = h(:);















