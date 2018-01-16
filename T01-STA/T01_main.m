% created on 2017-12-19
% Goal: do STA analysis for various cells
close all;
clc;clear;
set(0,'DefaultFigureWindowStyle','docked');

base_dir = 'C:\RathbumLab';

exp_dict =  T01_datalist();
for exp_id = exp_dict.keys()
    exp_id = char(exp_id);
    exp_data_dir = fullfile(base_dir,'Data\',exp_id,'\');
    for cell_id = exp_dict(exp_id)
        cell_id = char(cell_id);
        work_dir = fullfile(base_dir,'results\T01\',exp_id,'\',cell_id,'\');
        config_file = fullfile(exp_data_dir,'analysis_config.ini');

        if ~exist(work_dir,'dir'), mkdir(work_dir); end

        exp_ps = ini2struct(config_file); 

        exp_ps.exp_id = exp_id;
        exp_ps.cell_id = cell_id;
        exp_ps.work_dir = work_dir;
        exp_ps.data_dir = exp_data_dir;

        [STA_ps, D_ps] = STA_computation(exp_ps);
        T01_plots(STA_ps, D_ps, exp_ps);
        copyfile(config_file,work_dir);
        save(fullfile(work_dir,'workspace.mat'));
    end
end
