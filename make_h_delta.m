function h = make_h_delta(N, Nwindow)

%--------------------------------------------------------------------------------------------------
%
% fcxn:    make_h_delta.m
% descrip: constructs an ideal impulse with delay N (>=0) on interval [0: Nw-1]
% domain:  discrete
% inputs:  /N/        location of non-zero impulse
%          /Nwindow/  length of sampled impulse response
% output:  /h/        impulse response of ideal delay on [0: Nw-1]
% author:  JN Damask
%
%--------------------------------------------------------------------------------------------------

h = zeros(Nwindow, 1);
h(N) = 1;

