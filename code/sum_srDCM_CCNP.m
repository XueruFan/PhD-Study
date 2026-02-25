% =========================================================
% CCNP rDCM / srDCM summary script
% Windows version
% =========================================================

clear; clc;

% ===================== PATH =============================

data_root   = 'E:\PhDproject\Study3\output';
r_root      = fullfile(data_root, 'rDCM');
sr_root     = fullfile(data_root, 'srDCM');
stat_root   = fullfile(data_root, 'sum');

if ~exist(stat_root, 'dir')
    mkdir(stat_root);
end

nROI = 15;
nEC  = nROI * nROI;

%% ===================  rDCM  =============================

out_xlsx = fullfile(stat_root, 'CCNP_rDCM_summary.xlsx');

mat_files = dir(fullfile(r_root, '*_rDCM.mat'));
mat_files = mat_files(~startsWith({mat_files.name}, '._'));

nSub = numel(mat_files);
fprintf('Found %d CCNP rDCM files\n', nSub);

EC_all  = zeros(nSub, nEC);
subject = cell(nSub, 1);

for i = 1:nSub
    
    fn = mat_files(i).name;
    load(fullfile(r_root, fn));   % loads EC
    
    % 从文件名提取subject（去掉 _rDCM.mat）
    subject{i} = erase(fn, '_rDCM.mat');
    
    EC_all(i, :) = reshape(EC', 1, []);
    
end

% ===== 构造列名 =====
colNames = cell(1, nEC);
k = 1;
for from = 1:nROI
    for to = 1:nROI
        colNames{k} = sprintf('EC_%02d_to_%02d', from, to);
        k = k + 1;
    end
end

T = array2table(EC_all, 'VariableNames', colNames);
T = addvars(T, subject, 'Before', 1);

writetable(T, out_xlsx);
fprintf('CCNP rDCM summary Done!\n');


%% ===================  srDCM  ============================

out_xlsx = fullfile(stat_root, 'CCNP_srDCM_summary.xlsx');

mat_files = dir(fullfile(sr_root, '*_srDCM.mat'));
mat_files = mat_files(~startsWith({mat_files.name}, '._'));

nSub = numel(mat_files);
fprintf('Found %d CCNP srDCM files\n', nSub);

EC_all  = zeros(nSub, nEC);
subject = cell(nSub, 1);

for i = 1:nSub
    
    fn = mat_files(i).name;
    load(fullfile(sr_root, fn));   % loads EC_sparse
    
    subject{i} = erase(fn, '_srDCM.mat');
    
    EC_all(i, :) = reshape(EC_sparse', 1, []);
    
end

% ===== 构造列名 =====
colNames = cell(1, nEC);
k = 1;
for from = 1:nROI
    for to = 1:nROI
        colNames{k} = sprintf('EC_%02d_to_%02d', from, to);
        k = k + 1;
    end
end

T = array2table(EC_all, 'VariableNames', colNames);
T = addvars(T, subject, 'Before', 1);

writetable(T, out_xlsx);
fprintf('CCNP srDCM summary Done!\n');