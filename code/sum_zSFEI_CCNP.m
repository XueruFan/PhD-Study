clear; clc;

%% =========================
sfc_dir = 'E:\PhDproject\Study3\data\zSFEI';                 
out_dir = 'E:\PhDproject\Study3\output\zsfei';

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

nNet     = 15;
step_max = 7;

all_files = dir(fullfile(sfc_dir, '*_zSFEI_step*.mat'));
%%
for S = 1:step_max
    
    fprintf('Exporting step %d ...\n', S);
    
    step_tag   = sprintf('step%02d', S);
    step_files = all_files(contains({all_files.name}, step_tag));
    
    nSub = length(step_files);
    
    Site    = strings(nSub,1);
    Subject = strings(nSub,1);
    Session = strings(nSub,1);
    Run     = strings(nSub,1);
    Embed   = zeros(nSub, nNet);
    
    %% =========================
    for i = 1:nSub
        
        fname = step_files(i).name;
        fpath = fullfile(sfc_dir, fname);
        
        % 文件格式：
        % CCNPPEK0001_01_rest01_SFC_step01.mat
        
        tokens = regexp(fname, ...
            '^(CCNP[A-Z]+\d+)_(\d+)_(rest\d+)_DU15_zSFEI_step', ...
            'tokens');
        
        if isempty(tokens)
            warning('Filename not matched, skipped: %s', fname);
            continue
        end
        
        fullID  = tokens{1}{1};   % CCNPPEK0001
        Session(i) = tokens{1}{2};   % 01
        runStr  = tokens{1}{3};   % rest01
        
        %% ---- 从 fullID 中拆 Site 和 Subject ----
        id_tokens = regexp(fullID, ...
            '^CCNP([A-Z]+)(\d+)$', ...
            'tokens');
        
        Site(i)    = id_tokens{1}{1};   % PEK
        Subject(i) = id_tokens{1}{2};   % 0001
        
        Run(i) = erase(runStr, 'rest'); % 01
        
        load(fpath, 'zSFEI');
        Embed(i,:) = zSFEI(:)';
        
        clear zSFEI
    end

    varNames = ["Site", "Subject", "Session", "Run", ...
        arrayfun(@(x) sprintf("Net%02d", x), 1:nNet, 'UniformOutput', false)];
    
    T = array2table([Site, Subject, Session, Run, num2cell(Embed)], ...
        'VariableNames', varNames);
    
    out_name = sprintf('step%02d.xlsx', S);
    writetable(T, fullfile(out_dir, out_name));
    
end

fprintf('Done!\n');
