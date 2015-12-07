function target = make_sdt_target_far_lookback(bessel_order, length_pos_arm, length_neg_arm)

% function:   make_sdt_target_far_lookback
% descrip:    makes the target for 'far lookback'.
%
% inputs:     /bessel_order/        order of the bessel filters, taken as same degree below
%             /length_pos_arm/      M1 length on the + arm
%             /length_neg_arm/      M1 length on the - arm
%
%

bessel_def_pos = make_bessel_filter(bessel_order, length_pos_arm);
bessel_def_neg = make_bessel_filter(bessel_order, length_neg_arm);

% make the negative (longer) arm series and find the distant point at which the filter goes negative.
n_cand = [0: 10 * length_neg_arm];
neg_cand = plot_bessel_filter(bessel_def_neg, n_cand);
n_neg_crossover = find(neg_cand<=0, 2); 
N_filter = n_neg_crossover(2)-1;

% construct the pos and neg arms
s_pos = plot_bessel_filter(bessel_def_pos, [0:N_filter-1]);
s_neg = plot_bessel_filter(bessel_def_neg, [0:N_filter-1]);

S_pos = sum(s_pos);
S_neg = sum(s_neg);

% make correctly scaled impulse function, find index for zero cross-over
%   independently scale +/- components
h_target = s_pos(:) / S_pos - s_neg(:) / S_neg;
h_target(1) = 0;

%   scale combined +/- side
Sh_pos =  sum(h_target(h_target>0));
Sh_neg = -sum(h_target(h_target<0));
h_target = h_target / (mean([Sh_pos Sh_neg]));

% find zero crossing between pos and neg arms
N_crossover_pt = find(h_target(2:end-1)<0, 1, 'first') + 1;

% find pos / neg M1's
h_pos  = h_target(h_target>0);
h_neg  = h_target(h_target<0);
M1_pos = round([0:length(h_pos)-1] * h_pos);
M1_neg = length(h_pos) + round([0:length(h_neg)-1] * (-h_neg));

return_per_event_scale = 1 / (M1_neg - M1_pos);

% make delta-function target equivalent
h_delta_target = zeros(length(h_target), 1);
h_delta_target(M1_pos) =  1;
h_delta_target(M1_neg) = -1;

% make simple differencer with no lead or lag
h_cntl_target = zeros(length(h_target), 1);
h_cntl_target(1) = 1;
h_cntl_target(M1_neg-M1_pos) = -1;

% scale all impulse functions
h_target       = h_target       * return_per_event_scale;
h_delta_target = h_delta_target * return_per_event_scale;
h_cntl_target  = h_cntl_target  * return_per_event_scale;

% pack and return
target.h_target       = h_target;
target.h_delta_target = h_delta_target;
target.h_cntl_target  = h_cntl_target;
target.N_filter       = N_filter;
target.M1_pos         = M1_pos;
target.M1_neg         = M1_neg;
target.N_crossover_pt = N_crossover_pt;
target.return_per_event_scale = return_per_event_scale;



