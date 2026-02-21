clear; clc;

sfc_dir = '/Volumes/Zuolab_XRF/data/abide/zSFEI';       
out_dir = '/Volumes/Zuolab_XRF/output/abide/sfc/zsfei';

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

nNet     = 15;
step_max = 7;

%% =========================

all_files = dir(fullfile(sfc_dir, '*_zSFEI_step*.mat'));
all_files = all_files(~startsWith({all_files.name}, '._'));

for S = 1:step_max
    
    fprintf('Exporting step %d ...\n', S);
    
    step_tag   = sprintf('step%02d', S);
    step_files = all_files(contains({all_files.name}, step_tag));
    
    nSub = length(step_files);
    
    % 预分配
    Subject   = strings(nSub,1);
    Embed = zeros(nSub, nNet);
    
    % =========================
    for i = 1:nSub
        
        fname = step_files(i).name;
        fpath = fullfile(sfc_dir, fname);
        
        tokens = regexp(fname, '^(\d+)_', 'tokens');
        Subject(i) = tokens{1}{1};
        
        load(fpath, 'zSFEI');
        
        Embed(i,:) = zSFEI(:)';
        
        clear zSFEI
    end

    varNames = ["Subject", ...
        arrayfun(@(x) sprintf("Net%02d", x), 1:nNet, 'UniformOutput', false)];
    
    T = array2table([Subject, num2cell(Embed)], 'VariableNames', varNames);
    
    out_name = sprintf('step%02d.xlsx', S);
    writetable(T, fullfile(out_dir, out_name));
    
end

fprintf('Done!\n');
