function I12 = calc_overlap_integral(s1, s2)

%--------------------------------------------------------------------------------------------------
%
% fcxn:    calc_overlap_integral.m
% descrip: Computes the overlap integral between s1, s2 given that they are the
%          same length.
% inputs:  /s1/       a series vector
%          /s1/       a series vector
% output:  /I12/      normalized overlap integral
% author:  JN Damask
%
%--------------------------------------------------------------------------------------------------

% sanity
if length(s1) ~= length(s2),
    disp('s1 and s2 don''t have the the same lengths');
    I = -1;
    return;
end

% force
S1 = s1(:); S2 = s2(:);

% sums
I11 = sum(S1.^2);
I22 = sum(S2.^2);

% overlap integral
I12 = sum(S1.*S2) / sqrt(I11 * I22);








