% % created on 2017-12-19
% % Goal: do STA analysis for various cells
close all;
clc;clear;
set(0,'DefaultFigureWindowStyle','docked');

base_dir = 'D:\RathbumLab';

%exp_dict('2015_02_19') = {'adch_62f','adch_62g','adch_62l','adch_63h','adch_71b','adch_73h','adch_82f','adch_82k','adch_83d','adch_83e','adch_83f','adch_83g','adch_83h'};

exp_id = '2015_02_19';
cell_id = 'adch_62f';

exp_data_dir = fullfile(base_dir,'Data\',exp_id,'\');
work_dir = fullfile(base_dir,'results\T02\',exp_id,'\',cell_id,'\');
config_file = fullfile(exp_data_dir,'analysis_config.ini');

if ~exist(work_dir,'dir'), mkdir(work_dir); end

exp_ps = ini2struct(config_file); 

exp_ps.exp_id = exp_id;
exp_ps.cell_id = cell_id;
exp_ps.work_dir = work_dir;
exp_ps.data_dir = exp_data_dir;

[STA_ps, D_ps] = STA_computation(exp_ps);

%get_staIdx = @(splinedSTA_Idx) 2+fix((splinedSTA_Idx/((0.001/exp_ps.stimFreq)))+1);
get_staIdx = @(splinedSTA_Idx) 2+fix(splinedSTA_Idx/(((STA_ps.STA_t(2)-STA_ps.STA_t(1))/(STA_ps.splinedSTA_t(2)-STA_ps.splinedSTA_t(1)))+1));

sta_d1_idx = get_staIdx(D_ps.D1_idx);
sta_d2_idx = get_staIdx(D_ps.D2_idx);

if ~isnan(D_ps.D2_cross_ids(1)) && D_ps.D2_issig
    crop_i_idx = get_staIdx(D_ps.D2_cross_ids(1));
elseif ~isnan(D_ps.D1_cross_ids(1))&& D_ps.D1_issig
    crop_i_idx = get_staIdx(D_ps.D1_cross_ids(1));
else
    crop_i_idx = 1;
    display('Warning! No significant D2/D1 or no crossing were found! using the initial point of the STA.');
end

crop_j_idx = fix(length(STA_ps.STA)/2)+1;%ToDo - without adding one, the mid trial doesnt correspond to time zero

%% STA PLots
plot(STA_ps.STA_t,STA_ps.correctedSTA);hold on;
plot(STA_ps.STA_t(1,sta_d1_idx),STA_ps.correctedSTA(1,sta_d1_idx),'r*');
plot(STA_ps.STA_t(1,sta_d2_idx),STA_ps.correctedSTA(1,sta_d2_idx),'r*');

plot(STA_ps.STA_t(crop_i_idx:crop_j_idx),STA_ps.correctedSTA(crop_i_idx:crop_j_idx),'r--');

plt_ylim = [-1200, -400];
yaxis_line = zeros(length(plt_ylim(1):100:plt_ylim(2)));

plot(yaxis_line, plt_ylim(1):100:plt_ylim(2), 'k');

ylim([plt_ylim(1) plt_ylim(2)]);

STA_crop = STA_ps.correctedSTA(crop_i_idx:crop_j_idx);
STA_t_crop = STA_ps.STA_t(crop_i_idx:crop_j_idx);

for trialIdx = 1:size(STA_ps.tData,2)
    estim = STA_ps.tData(trialIdx).estim_amps';
    estim_ts = STA_ps.tData(trialIdx).estim_ts;
    estim_spts = STA_ps.tData(trialIdx).estim_spts;
    
    genSig = median(diff(STA_t_crop))*conv(STA_crop,estim,'full');
    genSig = genSig(1:length(estim_ts));
    
    figure(100+trialIdx);
    subplot(211);plot(estim_ts, estim);
    subplot(212);plot(estim_ts, genSig);
    xlim([1,100]);
    
    fig_basename = sprintf('%s_[%s]_trial=%d',exp_ps.exp_id,exp_ps.cell_id,trialIdx);
    subplot(211);title(fig_basename, 'Interpreter', 'none');
    saveas(gcf, [exp_ps.work_dir, fig_basename, '_STA_Conv2.jpeg']);
    break
end

%% Plot the histograms
figure;
estim_binedges = min(estim):50:max(estim);
[estim_bincounts] = histc(estim, estim_binedges);
subplot(121);bar(estim_binedges, estim_bincounts,'histc');

genSig_norm = (genSig-mean(genSig))/std(genSig);
genSig_binedges = min(genSig_norm):0.1:max(genSig_norm);
[genSig_bincounts] = histc(genSig_norm, genSig_binedges);
subplot(122);bar(genSig_binedges, genSig_bincounts,'histc');

%% spike trigering stimuli
% get the spike timings and while keeping certain values of estim and
% genSig before the spike time set all the other values to nan
% then plot the above figures again
figure(200+trialIdx);
subplot(211);plot(estim_ts, estim);

%subplot(212);plot(estim_t, genSig);
xlim([1,100]);
    
%T01_plots(STA_ps, D_ps, exp_ps);
%copyfile(config_file,work_dir)
