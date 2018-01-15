% % created on 2017-12-19
% % Goal: do STA analysis for various cells
close all;
clc;clear;
set(0,'DefaultFigureWindowStyle','docked');

base_dir = 'C:\RathbumLab';

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

crop_j_idx = fix(length(STA_ps.STA)/2);%the mid point can happen to not cross the exact zero point and that would be because we dont have samples there

%% Extracting STA and perparing variables
estim_meanline = STA_ps.estim_mean * ones(2 * exp_ps.tKerLen, 1);
STA_crop = STA_ps.correctedSTA(crop_i_idx:crop_j_idx);
STA_t_crop = STA_ps.STA_t(crop_i_idx:crop_j_idx);

Kw = length(STA_crop);

trialIdx = 15;

estim_amps = STA_ps.tData(trialIdx).estim_amps;
estim_inds = 1:length(estim_amps); % stimulus sample indices used for plotting
estim_ts = STA_ps.tData(trialIdx).estim_ts;
estim_spts = STA_ps.tData(trialIdx).estim_spts;
%% Figure 1xx - STA PLots
figsSeries = 100;
figure(figsSeries);

plot(STA_ps.STA_t,STA_ps.correctedSTA,'LineWidth',2);hold on;
plot(STA_ps.STA_t(1,sta_d1_idx),STA_ps.correctedSTA(1,sta_d1_idx),'r*');
plot(STA_ps.STA_t(1,sta_d2_idx),STA_ps.correctedSTA(1,sta_d2_idx),'r*');

plot(STA_t_crop,STA_crop,'r-','LineWidth',2);

plt_ylim = [-1200, -400];
yaxis_line = zeros(length(plt_ylim(1):100:plt_ylim(2)));

plot(yaxis_line, plt_ylim(1):100:plt_ylim(2), 'k');
plot(STA_ps.STA_t, estim_meanline, 'k');

ylim([plt_ylim(1) plt_ylim(2)]);
title(sprintf('STA plot for %s_[%s]',exp_ps.exp_id,exp_ps.cell_id), 'Interpreter', 'none');
%% Figure 1xx - XCorrelation of STA_crop with itset
figure(figsSeries+1);
sta_sta_xcorr_full = (1/exp_ps.stimFreq)*xcorr(STA_crop);

subplot(211);plot(1:Kw, STA_crop,'b');
subplot(212);plot(1:length(sta_sta_xcorr_full), sta_sta_xcorr_full,'b');
suptitle('Matlab Full XCorr of STA Crop with Itself');

figure(figsSeries+2);
sta_sta_xcorr_full2 = (1/exp_ps.stimFreq)*custom_xcorr(STA_crop,STA_crop,'full');

subplot(211);plot(1:Kw, STA_crop,'b');
subplot(212);plot(1:length(sta_sta_xcorr_full2), sta_sta_xcorr_full2,'b');
suptitle('Custom Full XCorr of STA Crop with Itself');

%% Figure 2xx - The Stimulus and the Generator Signal
figsSeries = 200;
figure(figsSeries+trialIdx);

genSig_vals = (1/exp_ps.stimFreq)*custom_xcorr(estim_amps,STA_crop);
genSig_inds = Kw:length(genSig_vals)+Kw-1;
genSig_ts = estim_ts(genSig_inds); %  We assign the timestamp corresponding to the end point of the xcorrel window to that genSig value

ax1 = subplot(211);plot(estim_inds, estim_amps);
ax2 = subplot(212);plot(genSig_inds, genSig_vals);
linkaxes([ax1,ax2],'x')
xlim([1,length(estim_amps)]);

fig_basename = sprintf('%s_[%s]_trial=%d',exp_ps.exp_id,exp_ps.cell_id,trialIdx);
subplot(211);title(fig_basename, 'Interpreter', 'none');
% saveas(gcf, [exp_ps.work_dir, fig_basename, '_STA_Conv2.jpeg']);

% plot(estim_ts, estim_amps,'LineWidth',4); hold on;
% plot(genSig_t, genSig, 'r','LineWidth',1);

%% Figure 3xx - Histograms
figsSeries = 300;
figure(figsSeries+trialIdx);

estim_binedges = linspace(min(estim_amps),max(estim_amps),50);
[estim_bincounts] = histc(estim_amps, estim_binedges);
subplot(121);bar(estim_binedges, estim_bincounts,'histc');

genSig_binedges = linspace(min(genSig_vals),max(genSig_vals),50);
[genSig_bincounts] = histc(genSig_vals, genSig_binedges);
subplot(122);bar(genSig_binedges, genSig_bincounts,'histc');

suptitle('Histogram of the Stimuli / Generator Signal Amplitudes')

%% Figure 4xx - Extracting Spike associated stimuli
% Extract the values in the stimuli and the generator signal that caused a
% spike. This would be a window of stimuli that immediately precede a spike
% or the single generator signal value before that spike
figsSeries = 400;
figure(figsSeries+trialIdx);

speriod = 1/exp_ps.stimFreq;%sampling period
pre_spike_sample = 16;% for 25 Hz, samples are 0.04 s far.

spike_estim_vals = interp1(estim_ts,estim_amps,estim_spts); % the values of the stimulus at the spike timepoint used for plotting
sp_assoc_stimuli = NaN(size(estim_amps));%spike associated stimuli
for spike_t = estim_spts' 
    idx_tochange = ((estim_ts>=(spike_t-pre_spike_sample*speriod))&(estim_ts<spike_t));
    sp_assoc_stimuli(idx_tochange) = estim_amps(idx_tochange);
end

ax1 = subplot(211);plot(estim_ts, estim_amps,'y');hold on;
plot(estim_ts, estim_amps,'k.');hold on;
plot(estim_spts, spike_estim_vals,'r.');hold on;
plot(estim_ts, sp_assoc_stimuli,'b*');hold on;

xlim([0,100]);

% In the first subplot the yellow line shows the stimuli variations,
% the black dots mark the sample points in the stimuli (might be overlayed by blue stars)
% the red dots show the time point that a spike occured
% the blue-stared values are the ones included as the spike associated stimuli
% note that in the variable sp_assoc_stimuli all the indices corresponding to black dots are NaN and indices correponding to blue stars are genSig

spike_genSig_vals = interp1(genSig_ts,genSig_vals,estim_spts);
sp_assoc_genSig = NaN(size(genSig_vals));%spike associated generator signal
for spike_t = estim_spts' 
    idx_tochange = ((genSig_ts>=(spike_t-speriod))&(genSig_ts<spike_t));
    sp_assoc_genSig(idx_tochange) = genSig_vals(idx_tochange);
end

ax2 = subplot(212);plot(genSig_ts, genSig_vals,'y');hold on;
plot(genSig_ts, genSig_vals,'k.');hold on;
plot(estim_spts, spike_genSig_vals,'r.');hold on;
plot(genSig_ts, sp_assoc_genSig,'b*');hold on;

xlim([0,100]);
linkaxes([ax1,ax2],'x')
suptitle('Extracting the spike associated stimuli/generator signal')

%% Figure 5xx - Histograms for the spike associated stimuli/generator signal
figsSeries = 500;
figure(figsSeries+trialIdx);

subplot(221);bar(estim_binedges, estim_bincounts,'histc');title('Stimuli');
subplot(222);bar(genSig_binedges, genSig_bincounts,'histc');title('Generator Signal');

nan_idx = isnan(sp_assoc_genSig);
sp_assoc_genSig_nonan = sp_assoc_genSig;
sp_assoc_genSig_nonan(nan_idx) = [];

estim_binedges = linspace(min(sp_assoc_stimuli),max(sp_assoc_stimuli),50);
[estim_bincounts] = histc(sp_assoc_stimuli, estim_binedges);
subplot(223);bar(estim_binedges, estim_bincounts,'histc');title('Spike Associated Stimuli');

sp_assoc_genSig_binedges = linspace(min(sp_assoc_genSig_nonan),max(sp_assoc_genSig_nonan),50);
[sp_assoc_genSig_bincounts] = histc(sp_assoc_genSig_nonan, sp_assoc_genSig_binedges);
subplot(224);bar(sp_assoc_genSig_binedges, sp_assoc_genSig_bincounts,'histc');title('Spike Associated Gnerator Signal');

suptitle('Histogram of the Spike Associated Stimuli / Generator Signal Amplitudes')

%% Figure 6xx - Firing Rate vs Generator Signal
% We would like to count the number of spikes corresponding to each value of the
% generator signal. For this we first assign a generator signal value to
% each spike time stamp. We then bin those spike_genSig_vals and count the
% number of times a value falls within each bin and visualize it as a histogram.
figsSeries = 600;
figure(figsSeries+trialIdx);

spike_genSig_vals = interp1(genSig_ts,genSig_vals,estim_spts(estim_spts>=genSig_ts(1)));
FRgenSig_binedges = linspace(min(spike_genSig_vals),max(spike_genSig_vals),50);
[FRgenSig_bincounts] = histc(spike_genSig_vals, FRgenSig_binedges);
bar(FRgenSig_binedges,FRgenSig_bincounts,'histc');
xlabel('Gen. Sig');
ylabel('#Spikes');

suptitle('Histogram of the #Spikes per bin of the Generator Signal Amplitude');
