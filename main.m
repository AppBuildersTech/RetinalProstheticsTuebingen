% created on 2017-12-19
% Goal: do STA analysis for various cells
close all;
clc;clear;
set(0,'DefaultFigureWindowStyle','docked');

exp_id = '2015_02_19';
cell_id = 'adch_62f';

data_dir = fullfile('D:\RathbumLab\Data\',exp_id,'\');
work_dir = fullfile('.\results\',exp_id,'\');
config_file = fullfile(data_dir,'analysis_config.ini');

if ~exist(work_dir,'dir'), mkdir(work_dir); end

exp_ps = ini2struct(config_file); 

exp_ps.exp_id = exp_id;
exp_ps.cell_id = cell_id;
exp_ps.work_dir = work_dir;
exp_ps.data_dir = data_dir;

[STA_ps, D_ps] = STA_computation(exp_ps);
STA_plots(STA_ps, D_ps, exp_ps);