function h = make_h_polyema_ct(tau, order, T, dt)

%--------------------------------------------------------------------------------------------------
%
% fcxn:    make_h_polyema_ct.m
% descrip: constructs a continuous-time poly-ema impulse response 
%          parameterized by /tau/ and sampled on the interval [0: dt: T-dt].
% domain:  continuous 
% note:    These ct poly-emas are the inverse Laplace xform of
%            H(s) = 1 / (1 + tau s)^order
%          which is to say, M1 = order x tau.
% inputs:  /tau/      tau, the first moment
%          /order/    order of pole
%          /T/        final time
%          /dt/       time interval
% output:  /h/        impulse response on [0:dt:T-dt] of poly-ema
% author:  JN Damask
%
%--------------------------------------------------------------------------------------------------

% sanity
if (order < 1) || (order > 5),
    disp('order out of bounds:   1 <= order <= 5');
    return
end

% dictionary of impulse responses,
hdict{1} = @(t)(1/tau^order * exp(-t/tau));
hdict{2} = @(t)(1/tau^order * t .* exp(-t/tau));
hdict{3} = @(t)(1/tau^order / 2 * t.^2 .* exp(-t/tau));
hdict{4} = @(t)(1/tau^order / 6 * t.^3 .* exp(-t/tau));
hdict{5} = @(t)(1/tau^order / 24 * t.^4 .* exp(-t/tau));

% make t axis
t = [0: dt: T-dt];

% compute
h = hdict{order}(t);
h = h(:);









