% =========================================================
% rDCM pipeline for fslr10k dataset
% =========================================================

clear; clc;

addpath(genpath('C:\matlab-toolbox\tapas'));

% ===================== PATH =============================

data_root = 'E:\PhDproject\Study3\data\timeseries\fslr10k';

out_root_r  = 'E:\PhDproject\Study3\output\rDCM';
out_root_sr = 'E:\PhDproject\Study3\output\srDCM';

if ~exist(out_root_r, 'dir'),  mkdir(out_root_r);  end
if ~exist(out_root_sr, 'dir'), mkdir(out_root_sr); end

% ===================== SETTINGS =========================

TR   = 2.0;
nROI = 15;

% ========================================================

ts_files = dir(fullfile(data_root, '*_network_ts.mat'));
fprintf('Found %d time series files\n', numel(ts_files));

for f = 1:numel(ts_files)

    fn = ts_files(f).name;
    full_fn = fullfile(data_root, fn);

    fprintf('\nProcessing: %s\n', fn);

    load(full_fn); 

    if size(ROI_ts,1) ~= nROI
        error('ROI number mismatch in %s', fn);
    end

    % -------- 保留完整原始前缀 --------
    % CCNPPEK0003_01_rest02_DU15_network_ts.mat
    % -> CCNPPEK0003_01_rest02_DU15

    baseName = erase(fn, '_network_ts.mat');

    % -------- 构造 Y --------

    Y.y  = ROI_ts';   % [T x ROI]
    Y.dt = TR;

    DCM = tapas_rdcm_model_specification(Y, [], []);

    % =====================================================
    % Classic rDCM
    % =====================================================

    rDCM = tapas_rdcm_estimate(DCM, 'r', [], 1);
    EC   = rDCM.Ep.A;

    save(fullfile(out_root_r, ...
        sprintf('%s_rDCM.mat', baseName)), ...
        'EC','TR');

    % =====================================================
    % Sparse rDCM
    % =====================================================

    rDCM_sparse = tapas_rdcm_estimate(DCM, 'r', [], 2);
    EC_sparse   = rDCM_sparse.Ep.A;

    save(fullfile(out_root_sr, ...
        sprintf('%s_srDCM.mat', baseName)), ...
        'EC_sparse','TR');

    % -------- 清理变量 --------
    clear ROI_ts_bilat Y DCM rDCM rDCM_sparse EC EC_sparse

end

fprintf('\nAll rDCM estimation finished!\n');