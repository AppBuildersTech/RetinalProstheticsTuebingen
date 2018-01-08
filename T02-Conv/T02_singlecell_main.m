% created on 2017-12-19
% Goal: do STA analysis for various cells
close all;
clc;clear;
set(0,'DefaultFigureWindowStyle','docked');

base_dir = 'C:\RathbumLab';

%exp_dict('2015_02_19') = {'adch_62f','adch_62g','adch_62l','adch_63h','adch_71b','adch_73h','adch_82f','adch_82k','adch_83d','adch_83e','adch_83f','adch_83g','adch_83h'};

exp_id = '2015_02_19';
exp_data_dir = fullfile(base_dir,'Data\',exp_id,'\');

cell_id = 'adch_62f';

work_dir = fullfile(base_dir,'results\',exp_id,'\',cell_id,'\');
config_file = fullfile(exp_data_dir,'analysis_config.ini');

if ~exist(work_dir,'dir'), mkdir(work_dir); end

exp_ps = ini2struct(config_file); 

exp_ps.exp_id = exp_id;
exp_ps.cell_id = cell_id;
exp_ps.work_dir = work_dir;
exp_ps.data_dir = exp_data_dir;

[STA_ps, D_ps] = STA_computation(exp_ps);
STA_plots(STA_ps, D_ps, exp_ps);
copyfile(config_file,work_dir)
