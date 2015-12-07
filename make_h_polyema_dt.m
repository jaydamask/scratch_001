function h = make_h_polyema_dt(M1_ref, order, Nwindow)

%--------------------------------------------------------------------------------------------------
%
% fcxn:    make_h_polyema_dt.m
% descrip: constructs a digital poly-ema impulse response 
%          parameterized by /M1_ref/ and /order/ on the interval [0: Nwindow-1].
% domain:  discrete
% note:    These dt poly-emas are the inverse z xform of
%            H(z) = (1-p)^order / (1 - pz/inverse)^order
%          which is to say, M1 = order x M1_ref.
% inputs:  /M1_ref/   first moment M1, for reference
%          /order/    order of pole
%          /Nwindow/  window length
%          /dt/       time interval
% output:  /h/        impulse response on [0: Nwindow-1] of poly-ema
% author:  JN Damask
%
%--------------------------------------------------------------------------------------------------

% sanity
if (order < 1) || (order > 5),
    disp('order out of bounds:   1 <= order <= 5');
    return
end

% p parameter, from M1_ref
p = M1_ref / (M1_ref + 1);

% gain adj
ginv = (1 - p)^order;

% dictionary of impulse responses
hdict{1} = @(n)( ginv * p.^n );
hdict{2} = @(n)( (n+1) .* hdict{1}(n) );
hdict{3} = @(n)( (n+2).*(n+1) / 2 .* hdict{1}(n) );
hdict{4} = @(n)( (n+3).*(n+2).*(n+1) / 6 .* hdict{1}(n) );
hdict{5} = @(n)( (n+4).*(n+3).*(n+2).*(n+1) / 24 .* hdict{1}(n) );

% make n axis
naxis = [0: Nwindow - 1];

% compute
h = hdict{order}(naxis);
h = h(:);











