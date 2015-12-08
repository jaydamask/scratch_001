function result = build_sdt_forecast_target(...
    strat_id, code, ref_date, target_scale, model_path, warehouse_path ...
    )

% function:   build_sdt_forecast_target
% descrip:
%
% inputs:     /strat_id/            the associated strat_id
%             /code/                instrument code
%             /ref_date/            YYYY.MM.DD format, a string of the ref_date for this data.
%             /target_scale/        a coefficient to the target scale, e.g. 2 (for 2x)
%             /model_path/          top-level path under which split EQ*.indicators.csv file(s) reside
%             /warehouse_path/      path to the data_warehouse
%
%  The input file is $model_path/signals/<strat_id>/<ref_date>/<code>.indicators.csv
%    and is not a q-converted file from indicators, but a direct split from the siggen code.
%
% warehouse:  /result_file/         table of target values
%
%  The output file is $warehouse_path/targets/<strat_id>/<code>.<ref_date>.targets.csv
%

% trace
result.trace = 'build_sdt_forecast_target';
mfn = mfilename;

% globals
global LOG_DEBUG;
global LOG_INFO;
global LOG_NOTICE;
global LOG_CRITICAL;

LOG_DEBUG=0;
LOG_INFO=1;
LOG_NOTICE=2;
LOG_CRITICAL=3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data, model parameters, inspect output mode

% EQ*.indicators.csv
%     timestamp,event,screen,rCnt,trd,sosc_0,sosc_1,sosc_2,mid,xlb,iym,spy,fxe
%     34257416785,81,39,1,0,0.5,0.5,0.5,21418.3,13735.9,20651.7,26311.8,27758.1
%     34257416785,81,39,2,0,0.5,0.5,0.5,21418.3,13734.4,20649.5,26311.4,27758.1
%     34257416785,81,39,3,0,0.5,0.5,0.5,21418.3,13734.4,20649.5,26311.8,27758.1
%     ..
%
% where I take rCnt and mid. These mids are 1e4 * log(px). 

% param_data:
%     name,val
%     target_0,25.0     <-- known to be smallest scale
%     target_1,57.4456
%     target_2,181.659
%     target_3,574.456

data_path   = [model_path '/signals/' strat_id '/' ref_date];
result_path = [warehouse_path '/targets/' strat_id '/' ref_date];
result_path = ['~/workspace/local_warehouse/sdt/targets/' strat_id '/' ref_date];
data_file   = [code '.indicators.csv'];
param_file  = [code '.target_length_scales.csv'];
result_file = [code '.targets.csv'];

logline(mfn, ['Data file   : ' data_path '/' data_file]);
logline(mfn, ['Params file : ' data_path '/' param_file]);
logline(mfn, ['Result file : ' result_path '/' result_file]);

model_data = importdata([data_path '/' data_file], ',', 1);
param_data = importdata([data_path '/' param_file], ',', 1);

% extract [rCnt,mid]
i_rCnt = strmatch('rCnt', model_data.colheaders);
i_mid  = strmatch('mid', model_data.colheaders);

% extract the mid and make control
log_mid  = model_data.data(:,i_mid);
sig_cntl = ones(size(log_mid));

% extract parameters
i_target = 0;
f = @(x,i)(strcmp(x, ['target_', num2str(i)]));
g = @(x)(f(x,i_target));
indx = find(cellfun(g, param_data.textdata(2:end,1)) == 1);
while ~isempty(indx),
   
    target_v(1+i_target) = round( ...
        param_data.data(find(cellfun(g, param_data.textdata(2:end,1)) == 1)) ...
        );
    
    i_target= i_target + 1;
    g = @(x)(f(x, i_target));
    indx = find(cellfun(g, param_data.textdata(2:end,1)) == 1);
    
end

% bessel M1 / Nx scaling
gamma_bessel = 1./ [1 3.63 2.7 2.29 2.07];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build target impulse response across target_indices

FAR_LOOK_AHEAD = 0;
FAR_LOOK_BACK = 1;
SMOOTH_DOMAIN = 2;

BESSEL_ORDER = 3;
% TARGET_TYPE = FAR_LOOK_BACK;
TARGET_TYPE = SMOOTH_DOMAIN;

targ_m = zeros(length(log_mid), length(target_v)-1);
for i_target = 2 : length(target_v),

    switch TARGET_TYPE
        
        case FAR_LOOK_AHEAD
            
            it = i_target - 1;
            target_pos_arm = round( target_v(i_target) * target_scale );
            target_neg_arm = round( target_v(1) / 2 ); % factor of 2 attempts to account for the width of the impulse rather than its M1
            
            logline(mfn, ['building target(' num2str(it) ') with arm lengths (pos,neg): (' ...
                num2str(target_pos_arm), ',' ...
                num2str(target_neg_arm) ')']);
            
            bessel_def_pos = make_bessel_filter(BESSEL_ORDER, target_pos_arm);
            bessel_def_neg = make_bessel_filter(BESSEL_ORDER, target_neg_arm);
            
            % make the positive arm series and find the distant point at which the filter goes negative.
            n_cand = [0: 10 * target_pos_arm];
            s_cand = plot_bessel_filter(bessel_def_pos, n_cand);
            n_crossover = find(s_cand<=0, 2); N_max = n_crossover(2)-1;
            
            % construct the pos and neg arms
            s_pos = plot_bessel_filter(bessel_def_pos, [0:N_max-1]);
            s_neg = plot_bessel_filter(bessel_def_neg, [0:N_max-1]);
            
            Spos = sum(s_pos);
            Sneg = sum(s_neg);
            
            % make correctly scaled window function, find index for minimum
            w_target = s_pos(:) / Spos - s_neg(:) / Sneg;
            [tmp, N_neg_mode] = min(w_target);
            
            n_ruler = [0: N_max-1] - N_neg_mode + 1;
            
            % now, of course a window and impulse response are flips of one another
            h_target = flipud(w_target);
            
            % keep a control so I know exactly how to associate the target with the original mid series
            h_cntl = [zeros(N_max - N_neg_mode,1); 1];
            
            % build training targets -- exec convolutions
            cand      = conv(h_target, log_mid - log_mid(1));
            cand_cntl = conv(h_cntl, sig_cntl);
            
            % select from cand data that is aligned to original
            targ_m(:, it) = cand(cand_cntl>0);
            targ_m(end-(N_max-N_neg_mode)+1:end, it) = NaN;
            
            % cache
            h_target_cache{it} = h_target(~isnan(h_target));
            
        case FAR_LOOK_BACK
            
            it = i_target - 1;
            target_pos_arm = round( target_scale * min([target_v(1) target_v(i_target) / 20])  ); 
            target_neg_arm = round( target_v(i_target) * target_scale );
            
            % make the impulse response
            T = make_sdt_target_far_lookback(BESSEL_ORDER, target_pos_arm, target_neg_arm);

            logline(mfn, ['built target(' num2str(it) ') with arm M1 lengths (pos,neg): (' ...
                num2str(T.M1_pos), ',' ...
                num2str(T.M1_neg) ')']);
            
            % exec convolution
            cand       = conv(T.h_target, log_mid - log_mid(1));
            cand_delta = conv(T.h_delta_target, log_mid - log_mid(1));
            cand_cntl  = conv(T.h_cntl_target, log_mid - log_mid(1));
            
            % carefully account for non-causal component of target
            cand_trim       = cand(T.N_filter: length(log_mid));
            cand_delta_trim = cand_delta(T.N_filter: length(log_mid)); 
            cand_cntl_trim  = cand_cntl(T.N_filter: length(log_mid)); 

            t_vec = NaN * ones(length(log_mid),1);
            t_vec(T.N_filter - T.N_crossover_pt + 1: end-T.N_crossover_pt+1) = cand_trim;

            targ_m(:, it) = t_vec;
                        
            % cache
            Tc{it} = T;
            
        case SMOOTH_DOMAIN
            
            it = i_target - 1;
            
            lag_scale  = 1.0 - 0.5 * target_v(i_target) / target_v(end);
            ref_scale  = 1.0 - lag_scale;
            lead_scale = ref_scale * (1 - gamma_bessel(BESSEL_ORDER));
            
            target_pos_arm = round( lead_scale * target_v(end) ) * target_scale;
            target_neg_arm = round( lag_scale  * target_v(end) ) * target_scale;
            
            % make the impulse response
            T = make_sdt_smoothdomain_target(BESSEL_ORDER, target_pos_arm, target_neg_arm);
            
            logline(mfn, ['built target(' num2str(it) ') with arm ref positions (lead,lag): (' ...
                num2str(T.rtn_position_ref(1)), ',' ...
                num2str(T.rtn_position_ref(2)) ') ' ...
                'N_filter: ' num2str(diff(T.rtn_position_ref))]);
            
            % exec convolution
            cand       = conv(T.h_target, log_mid - log_mid(1));
            cand_delta = conv(T.h_delta_target, log_mid - log_mid(1));
            cand_cntl  = conv(T.h_cntl_target, log_mid - log_mid(1));
            
            % carefully account for non-causal component of target
            cand_trim       = cand(T.N_filter: length(log_mid));
            cand_delta_trim = cand_delta(T.N_filter: length(log_mid)); 
            cand_cntl_trim  = cand_cntl(T.N_filter: length(log_mid)); 

            t_vec = NaN * ones(length(log_mid),1);
            t_vec(T.N_filter - T.N_crossover_pt + 1: end-T.N_crossover_pt+1) = cand_trim;

            targ_m(:, it) = t_vec;
                        
            % cache
            Tc{it} = T;
                        
    end
    
end

% some logging
logline(mfn, ['series length: ' num2str(length(log_mid))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build dataset and export

header={'rCnt', 'log_mid'};
for k = 2: length(target_v),
    %header = cat(2, header, ['t_' num2str(target_v(k)) ]);
    header = cat(2, header, ['targ_' num2str(k-2) ]);
end

logline(mfn, ['exporting ' num2str(size(log_mid,1)) ' records']);
if (~isdir(result_path))
    if (~mkdir(result_path))
        disp(['ERROR -- cannot mkdir on path: ' result_path]);
        return
    end
end

matrix_limited_export(...
    [model_data.data(:,i_rCnt) log_mid(:) targ_m], ...
    header, ...
    [result_path '/' result_file] ...
    );

% for yucks
result.log_mid = log_mid;
result.targ_m  = targ_m;
result.targ_v  = target_v;
result.Tc      = Tc;





