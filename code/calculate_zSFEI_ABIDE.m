clear; clc;

%% =========================
data_root = '/Volumes/Zuolab_XRF/data/abide/timeseries';
out_root  = '/Volumes/Zuolab_XRF/data/abide/zSFEI';

if ~exist(out_root, 'dir')
    mkdir(out_root);
end

addpath(genpath('/Users/xuerufan/matlab-toolbox/CCS-master'));  % CCS 工具箱路径
addpath(genpath('/Users/xuerufan/PhD-Study3/code'));
cd('/Users/xuerufan/PhD-Study3/code')

files = dir(fullfile(data_root, '*_network_ts.mat'));
files = files(~startsWith({files.name}, '._'));

nSub  = length(files);

step_max  = 6;      % 最大步数

%% =========================
for iSub = 1:nSub
    
    fprintf('Subject %d / %d\n', iSub, nSub);
    
    load(fullfile(data_root, files(iSub).name));
    
    N = size(ROI_ts,1); % ROI个数
    
    % -------- FC --------
    FC = corr(ROI_ts');
    FC(1:N+1:end) = 0; %自相关为0
    
    absFC = abs(FC);
    
    % -------- 稀疏化：每个网络至少 1 条边 --------
    adj = zeros(N, N);
    
    for i = 1:N
        row = absFC(i, :);
        row(i) = 0;
        [~, idx] = max(row);
        adj(i, idx) = 1;
    end
    
    % 无向化
    adj = max(adj, adj');

    for S = 1:step_max
        % -------- SFC --------
        SFC = ccs_core_graphwalk_old(adj, S, 'nbtw');
        SFC(1:N+1:end) = 0;
    
        % -------- SFEI --------
        tmp = SFC;
        tmp(1:N+1:end) = NaN;
        SFEI = mean(tmp, 2, 'omitnan');
    
        % -------- zSFEI --------
        muS  = mean(SFEI);
        sigS = std(SFEI);
    
        zSFEI = (SFEI - muS) ./ sigS;
        zSFEI(isnan(zSFEI)) = 0;
    
        [~, oldname, ~] = fileparts(files(iSub).name);
        tokens = regexp(oldname, '^(\d+)_', 'tokens');
        subID  = tokens{1}{1};
            
        save_name = sprintf('%s_DU15_zSFEI_step%02d.mat', subID, S);
        save(fullfile(out_root, save_name), 'SFC', 'SFEI', 'zSFEI');
    
        clear SFC SFEI zSFEI muS sigS tmp
    end
    clear ts FC absFC adj
end

fprintf('Done!\n');
