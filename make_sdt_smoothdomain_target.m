function target = make_sdt_smoothdomain_target(bessel_order, M1_lead_arm, M1_lag_arm)


% function:   make_sdt_smoothdomain_target
% descrip:    makes a 'smooth-domain' target, which is to say, a target with a causal
%             and non-causal component, where the causal component smooths (the input, which
%             is (log)px) and the non-causal component seeks to forecast the motion of the
%             smooth.
%
%             The impulse response /h_targ/ looks like this:
%
%                          ^ (-rtn_pos_ref)
%                          |
%                 ++++++++++++++++++|
%             ----------------------0--------------------------------------------> n
%                                    ---------------------------------------|
%                                                        |
%                                                        v (+rtn_neg_ref)
%
%                 (<    N+_max    >)0(<               N-_max               >)
%                 |--    N_x     -->|
%                 |--                   N_filter                         -->|
%                 |--   -->|  M1_lead
%                                   |--               -->|  M1_lag
%
%
%             where the gain is zero, the +/- arm integrals are 1 / (M1_neg - M1_pos), 
%             and the +/- cut-offs are the natural zero crossing of the bessel filters.
%
% inputs:     /bessel_order/   order of the bessel filters, taken as same degree below
%             /M1_lead_arm/    M1 length on the lead, + arm
%             /M1_lag_arm/     M1 length on the lag,- arm
%
%

% get bessel definitions
bessel_def_lead = make_bessel_filter(bessel_order, M1_lead_arm);
bessel_def_lag  = make_bessel_filter(bessel_order, M1_lag_arm);

% make lead arm
n_cand = [0: 10 * M1_lead_arm];
cand = plot_bessel_filter(bessel_def_lead, n_cand);
cand(1) = max([cand(1) 0]);

N_crossover_lead = find(cand<0, 1, 'first');
h_lead = plot_bessel_filter(bessel_def_lead, [0: N_crossover_lead-1]);
h_lead = h_lead(:) / sum(h_lead); 

% make lag arm
n_cand = [0: 10 * M1_lag_arm];
cand = plot_bessel_filter(bessel_def_lag, n_cand);
cand(1) = max([cand(1) 0]);

N_crossover_lag = find(cand<0, 1, 'first');
h_lag = plot_bessel_filter(bessel_def_lag, [0: N_crossover_lag-1]);
h_lag = h_lag(:) / sum(h_lag); 

% compute return-per-event scale
M1_lead = round([0: length(h_lead)-1] * h_lead);
M1_lag  = round([0: length(h_lag)-1] * h_lag);

M1_lead_wrt_zero = N_crossover_lead - M1_lead;
M1_lag_wrt_zero  = M1_lag;

return_per_event_scale = 1 / (M1_lead_wrt_zero + M1_lag_wrt_zero);

% combine lead and lag into a single impulse response /h/
h_target = [h_lead; -h_lag];
N_x      = N_crossover_lead;
N_filter = length(h_target);

% make delta-function target equivalent
h_delta_target = zeros(N_filter, 1);
h_delta_target(M1_lead)      =  1;
h_delta_target(N_x + M1_lag) = -1;

% make simple differencer with no lead or lag
h_cntl_target = zeros(N_filter, 1);
h_cntl_target(1) = 1;
h_cntl_target(N_x + M1_lag - M1_lead) = -1;

% scale all responses
h_target       = h_target * return_per_event_scale;
h_delta_target = h_delta_target * return_per_event_scale;
h_cntl_target  = h_cntl_target  * return_per_event_scale;

% pack for return
target.h_target               = h_target;
target.h_delta_target         = h_delta_target;
target.h_cntl_target          = h_cntl_target;
target.N_filter               = N_filter;
target.N_crossover_pt         = N_x;
target.M1                     = [M1_lead M1_lag];
target.rtn_position_ref       = [-M1_lead_wrt_zero M1_lag_wrt_zero];
target.return_per_event_scale = return_per_event_scale;





