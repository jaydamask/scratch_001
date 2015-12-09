function result = process_panel_regression(...
    set_name, ref_date, code, target_name, redecimation_interval, warehouse_path ...
    )

% function:   process_panel_regression
% descrip:    Loads specified panel and performs OLS. Saves result back to d-warehouse.
% update:     The panel is not fully decimated. Complete decimation is where incoming panel
%             is decimated at the /redecimation_interval/. This code walks thru integral offsets
%             across the /redecimation_interval/ and aggregates the beta (etc) results.
%
% inputs:     /ref_date/              YYYYMMDD format, a string of the ref_date for this data.
%             /code/                  instrument code
%             /target_name/           name of the target as found on the file
%             /redecimation_interval/ the incoming panel is further decimated at this interval
%             /warehouse_path/        path to the data_warehouse
%
%  The input file is $path/<ref_date>.<code>.<target_name>.panel.csv
%
% warehouse:  /result_file/         report of results
%             /beta_file/           report of betas
%
%  The output file is $path/<ref_date>.<code>.<target_name>.regression.csv
%                     $path/<ref_date>.<code>.<target_name>.beta.csv
%

% trace
result.trace = 'process_panel_regression';
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

% signals data
%     s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32,s33,s34,s35,s36,s37,s38,s39,s40,s41,s42,s43,s44,s45,s46,s47,s48,s49,s50,s51,s52,s53,s54,s55,s56,s57,s58,s59,t0,t1,t2,t3,t4
%     0.187657,0.395618,0.221276,0.282603,0.28773,-1.17245,-0.698216,-0.558219,0.0214623,-0.000923835,0.828102,0.686969,0.635741,0.363614,0.412631,0.261811,-0.0271532,0.166789,1.1774,1.07625,-4.33864,0.365224,-2.93225,-9.21268,-8.96769,-5.01798,-7.3047,-9.72786,-17.2808,-17.1177,-1.73515e-08,-0.0162572,-0.00666837,0.0013797,-1.83267e-05,-0.0390658,-4.01785,-0.655975,0.0886558,-0.00214055,-0.341087,-0.096983,-0.0620613,-0.0131235,0.000148445,-0.341087,-0.095968,-0.0658578,-0.0123182,0.000132952,-0.341087,-0.0959717,-0.078521,-0.0107571,0.000101638,0,0,0,0,0,-0.000206835,-0.001203924,-0.002564516,-0.002102959,-0.002816459
%     ..

data_path    = [warehouse_path '/' set_name '/data'];
result_path  = [warehouse_path '/' set_name '/result'];
result_path  = ['~/workspace/local_warehouse/targets/result'];

data_file    = [code '.' ref_date '.' target_name '.panel.csv'];
beta_file    = [code '.' ref_date '.' target_name '.beta.csv'];
summary_file = [code '.' ref_date '.' target_name '.is.summary.csv'];
offset_file  = [code '.' ref_date '.' target_name '.offset.summary.csv'];

logline(mfn, ['Data file    : ' data_path '/' data_file]);
logline(mfn, ['Beta file    : ' result_path '/' beta_file]);
logline(mfn, ['Summary file : ' result_path '/' summary_file]);
logline(mfn, ['Offset file  : ' result_path '/' offset_file]);

panel = importdata([data_path '/' data_file],',',1);

% the following function won't work because there are fileds in between , are null
% /dat/hftanalytic/calculation-alpha/data_warehouse/splp/data/EQ0010014700001000.20110331.t0.panel.csv
% [panel.data panel.colheaders] = load_csvfile_one_header([data_path '/' data_file]);

% extract target column
i_targ = find(cellfun(@(x)(strcmp(x, target_name)), panel.colheaders) == 1);
if isempty(i_targ),
   logline(mfn, 'ERROR: cannot find the target column in the panel');
   return
else
   logline(mfn, ['target panel column: ' num2str(i_targ)]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build regression

% some logging
logline(mfn, ['panel length: ' num2str(length(panel.data))]);

% some hard-coding
Nspd = 4;
Npoly = 4;
Nqtind = 5;
signal_vector = 1: Nqtind * Npoly + Nspd * Npoly;
if (regexp(set_name,'XLEOIH','once') == 1)
    signal_vector = 1:30; % instead of 35;
elseif (regexp(set_name,'XLB','once') == 1)
     signal_vector = 1:36;
elseif (regexp(set_name,'OIH','once') == 1)   
    signal_vector = 1:30;
else
    disp(['invalid set_name' set_name]);
    return
end

% sigind_vector = [
%     1 % cntden 0
%     1 % cntden 1
%     1 % cntden 2
%     0 % bkszpx 0
%     0 % bkszpx 1
%     0 % bkszpx 2
%     1 % trd    1
%     0 % trd    2
%     1 % spd    1
%     1 % spd    2
%     1 % spd    3
%     1 % spd    4
% ];
% sind = repmat(sigind_vector, 1, 4)'; sind = sind(:);
% signal_vector = find(sind>0);

% make target and regressors

% iterate offsets thru the redecimation_interval
valid_vec = find(~isnan(panel.data(:, i_targ)));
valid_panel = panel.data(valid_vec, :);
for offset = 1: redecimation_interval,


    Panel = valid_panel(offset: redecimation_interval: end, :);
    targ  = Panel(:, i_targ);

%     Panel = panel.data(offset: redecimation_interval: end, :);
%     targ = Panel(:, i_targ);

    X = ones(size(Panel,1), 1 + length(signal_vector));
    X(:, 2:end) = Panel(:, signal_vector);
    
    [beta, betatol, R2, res, dw] = simple_regression(targ, X);

    beta_m(:,offset) = beta;
    R2_v(offset)     = R2;
    
    % Z regression
    Z = (X - ones(size(X)) * diag(mean(X))) * diag(1./std(X));
    Z(:,1) = 1;
    
    C = Z(:,2:end)' * Z(:,2:end) / size(Z,1);
    
    [zbeta, zbetatol, zR2, zres, zdw] = simple_regression(targ, Z);

    zbeta_m(:,offset) = zbeta;
    
    % some logging
    logline(mfn, ['R2: X - ' num2str(R2) '  Z - ' num2str(zR2)]);
    logline(mfn, ['dw: X - ' num2str(dw(1)) '  Z - ' num2str(zdw(1))]);
    
    % a plot (in interactive mode)
    figure(32); clf; grid on

    stem(zbeta + zbetatol, 'r'); hold on; grid on; stem(zbeta - zbetatol, 'r'); stem(zbeta, 'b')
    ylim([-1 1] * 2e-4); xlim([0 50]);
    title(['date:' ref_date ', offset: ' num2str(offset)]);
    drawnow
    
    figure(23); clf;
    surf(blkdiag(C,0));
    view(0,-90); axis equal; caxis([-0.25 1]);
    drawnow
        
end

% aggregate betas
beta  = median(beta_m, 2);
zbeta = median(zbeta_m, 2);
R2    = median(R2_v);

% store for return val
result.beta_m = beta_m;
result.zbeta_m = zbeta_m;
result.R2_v = R2_v;
result.beta    = beta;
result.zbeta   = zbeta;
result.R2      = R2;
result.betatol = betatol;
result.zbetatol = zbetatol;
result.signal_vector = signal_vector;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build dataset and export

ds = dataset( ...
    {beta, 'beta'}, ...
    {zbeta, 'zbeta'}, ...
    {betatol, 'betatol'}, ...
    {zbetatol, 'zbetatol'});

logline(mfn, ['exporting ' num2str(length(beta)) ' records']);
export(ds, 'file', [result_path '/' beta_file], 'delimiter', ',');

ds2 = dataset( ...
    {R2, 'R2'}, ...
    {dw(1), 'dw'}, ...
    {std(targ), 'std_targ'}, ...
    {std(res), 'std_res'});

logline(mfn, ['exporting summary records']);
export(ds2, 'file', [result_path '/' summary_file], 'delimiter', ',');

ds3 = dataset( ...
    {R2_v', 'R2'}, ...
    {zbeta_m', 'zbeta'});

logline(mfn, ['exporting offset records']);
export(ds3, 'file', [result_path '/' offset_file], 'delimiter', ',');








