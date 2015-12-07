function result = build_forecast_target(...
    ref_date, code, warehouse_path, params_file, ...
    set_name ...
    )
    
% function:   build_forecast_target
% descrip:    
%
% inputs:     /ref_date/            YYYYMMDD format, a string of the ref_date for this data.
%             /code/                instrument code
%             /warehouse_path/      path to the data_warehouse
%             /params_file/         file that contains the parameters for this calculation
%
%  The input file is $path/<ref_date>.<code>.mids.csv
%
% warehouse:  /result_file/         report of forecasted daily_intensity
%
%  The output file is $path/<model>.<set_name>.<ref_date>.<code>.di_forecasted.csv
%

% trace
result.trace = 'build_forecast_target';
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

% mid data
%     rcnt,mid
%     2,58.775
%     3,58.775
%     4,58.775
%     ..

% param_data:
%     param,val
%     neff_window_ratio,0.5
%     window1,500
%     window2,1000
%     window3,2000
%     window4,3000
%     window5,4000

data_path   = [warehouse_path '/' set_name '/data'];
result_path = [warehouse_path '/' set_name '/result'];
result_path = ['~/workspace/local_warehouse/targets' '/' set_name '/result'];
data_file   = [code '.' ref_date '.mids.csv'];
result_file = [code '.' ref_date '.targets.csv'];

logline(mfn, ['Data file   : ' data_path '/' data_file]);
logline(mfn, ['Params file : ' data_path '/' params_file]);
logline(mfn, ['Result file : ' result_path '/' result_file]);

model_data = importdata([data_path '/' data_file],',',1);
param_data = importdata([data_path '/' params_file],',',1);

% extract params
indx = find(cellfun(@(x)(strcmp(x, 'neff_window_ratio')), param_data.textdata(2:end,1)) == 1);
if ~isempty(indx)
   neff_ratio =  param_data.data(indx);
end

i_window = 1;
f = @(x,i)(strcmp(x, ['window', num2str(i)]));
g = @(x)(f(x,i_window));
indx = find(cellfun(g, param_data.textdata(2:end,1)) == 1);
while ~isempty(indx),
   
    window_v(i_window) = param_data.data(find(cellfun(g, param_data.textdata(2:end,1)) == 1));
    
    i_window = i_window + 1;
    g = @(x)(f(x, i_window));
    indx = find(cellfun(g, param_data.textdata(2:end,1)) == 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build targets

% some logging
logline(mfn, ['series length: ' num2str(length(model_data.data))]);

log_mid_v = log(model_data.data(:,2));
sig_cntl = ones(size(log_mid_v));

targ_m = zeros(length(log_mid_v), length(window_v));
for iw = 1: length(window_v),
    
    logline(mfn, ['  target length: ' num2str(window_v(iw))]);
    
    % define target impulse function
    n = [0: window_v(iw)-1];
    Neff = neff_ratio * window_v(iw);
    d = Neff / (1 + Neff);
    iir = (1 - d) * d.^n;
    iir = iir / sum(iir);
    
    h_targ = [iir -1];
    
    % convolve with log_mid_v
    cand = conv(h_targ, log_mid_v - log_mid_v(1));
    
    % keep a control so I know exactly how to associate the target with the original mid series
    h_cntl = [zeros(size(iir)) 1];   
    cand_cntl = conv(h_cntl, sig_cntl); 
    
    % select from cand aligned to original data, 
    %   tag the end of the target vector where the impulse function extends beyond the bound
    targ_m(:,iw) = cand(cand_cntl>0);
    targ_m(end-length(h_targ)+1:end, iw) = NaN;
       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build dataset and export

header={'rcnt', 'log_mid'};
for k = 1: length(window_v),
    header = cat(2, header, ['t' num2str(window_v(k))]);
end

logline(mfn, ['exporting ' num2str(size(log_mid_v,1)) ' records']);
matrix_limited_export([model_data.data(:,1) log_mid_v(:) targ_m], header, [result_path '/' result_file])






