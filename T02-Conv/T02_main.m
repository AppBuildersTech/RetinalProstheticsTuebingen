% % created on 2017-12-19
% % Goal: do STA analysis for various cells
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

        crop2_idx = fix(length(STA_ps.STA)/2);%the mid point can happen to not cross the exact zero point and that would be because we dont have samples there

        if ~isnan(D_ps.D2_cross_ids(1)) && D_ps.D2_issig
            crop1_idx = get_staIdx(D_ps.D2_cross_ids(1));
        elseif ~isnan(D_ps.D1_cross_ids(1))&& D_ps.D1_issig
            crop1_idx = get_staIdx(D_ps.D1_cross_ids(1));
        elseif ~isnan(D_ps.D2_finsig_ids(1)) && D_ps.D2_issig
            crop1_idx = get_staIdx(D_ps.D2_finsig_ids(1));
        elseif ~isnan(D_ps.D1_finsig_ids(1))&& D_ps.D1_issig
            crop1_idx = get_staIdx(D_ps.D1_finsig_ids(1));
        else
            crop1_idx = 1;
            display('Warning! No significant D2/D1 or no crossing were found! using the initial point of the STA.');
        end

        %% Extracting STA and perparing variables
        estim_meanline = STA_ps.estim_mean * ones(2 * exp_ps.tKerLen, 1);
        STA_crop = STA_ps.correctedSTA(crop1_idx:crop2_idx);
        STA_t_crop = STA_ps.STA_t(crop1_idx:crop2_idx);

        Kw = length(STA_crop);

        %% Figure 1xx - STA PLots
        figIdx = 1;
        figure();

        subplot(4,1,[1,2]);
        plot(STA_ps.STA_t,STA_ps.correctedSTA, 'b', 'LineWidth',2);hold on;
        plot(STA_ps.STA_t(1,sta_d1_idx),STA_ps.correctedSTA(1,sta_d1_idx),'r*');
        plot(STA_ps.STA_t(1,sta_d2_idx),STA_ps.correctedSTA(1,sta_d2_idx),'r*');

        plot(STA_t_crop,STA_crop,'r-','LineWidth',2);

        plt_ylim = [-1200, -400];
        yaxis_line = zeros(length(plt_ylim(1):100:plt_ylim(2)));

        plot(yaxis_line, plt_ylim(1):100:plt_ylim(2), 'k');
        plot(STA_ps.STA_t, estim_meanline, 'k');

        ylim([plt_ylim(1) plt_ylim(2)]);
        title(sprintf('STA plot for %s_[%s]',exp_ps.exp_id,exp_ps.cell_id), 'Interpreter', 'none');

        subplot(413);plot(1:Kw, STA_crop,'b'); title('STA Crop, a.k.a Kernel');

        sta_sta_xcorr_full = (1/exp_ps.stimFreq)*custom_xcorr(STA_crop,STA_crop,'full');

        subplot(414);plot(1:length(sta_sta_xcorr_full), sta_sta_xcorr_full,'k-');title('Padded XCorrelation of Kernel with itset');

        saveas(gcf, [exp_ps.work_dir, sprintf('%s_[%s]',exp_ps.exp_id,exp_ps.cell_id), sprintf('_F%.2d.jpeg',figIdx)]);

        %% Rest of the figures

        for trialIdx = 17%1:length(STA_ps.tData)
            fig_basename = sprintf('%s_[%s]_trialIdx=%.2d',exp_ps.exp_id,exp_ps.cell_id,trialIdx);

            estim_amps = STA_ps.tData(trialIdx).estim_amps;
            estim_inds = 1:length(estim_amps); % stimulus sample indices used for plotting
            estim_ts = STA_ps.tData(trialIdx).estim_ts;
            estim_spts = STA_ps.tData(trialIdx).estim_spts;

            %% Figure 2xx - The Stimulus and the Generator Signal
            figIdx = 2;
            figure();
            
            estim_amps_norm = (estim_amps - STA_ps.estim_mean) / STA_ps.estim_std;
            STA_crop_norm = (STA_crop - STA_ps.estim_mean) / STA_ps.estim_std;
            
            genSig_vals = (1/exp_ps.stimFreq)*custom_xcorr(estim_amps_norm,STA_crop_norm);
            genSig_inds = Kw:length(genSig_vals)+Kw-1;
            genSig_ts = estim_ts(genSig_inds); %  We assign the timestamp corresponding to the end point of the xcorrel window to that genSig value

            ax1 = subplot(211);plot(estim_inds, estim_amps_norm, 'b');
            ax2 = subplot(212);plot(genSig_inds, genSig_vals,'b');
            linkaxes([ax1,ax2],'x');
            xlim([1,length(estim_amps)]);

            figTitle = sprintf('%s [%s]\n trialIdx = %.2d - The normlzd Stmuli and the Gen Signal computed for normlzd Stim/Kernel',strrep(exp_ps.exp_id,'_','.'),strrep(exp_ps.cell_id,'_','-'), trialIdx);
            suptitle(figTitle);
            saveas(gcf, [exp_ps.work_dir, fig_basename, sprintf('_F%.2d.jpeg',figIdx)]);

            %% Figure 3xx - Extracting Spike associated stimuli
            % Extract the values in the stimuli and the generator signal that caused a
            % spike. This would be a window of stimuli that immediately precede a spike
            % or the single generator signal value before that spike
            figIdx = 3;
            figure();

            speriod = 1/exp_ps.stimFreq;%sampling period
            pre_spike_sample = 16;% for 25 Hz, samples are 0.04 s far.

            spike_estim_vals = interp1(estim_ts,estim_amps_norm,estim_spts); % the values of the stimulus at the spike timepoint used for plotting
            sp_assoc_stimuli = NaN(size(estim_amps_norm));%spike associated stimuli
            sp_assoc_stimuli_weighted = [];
            for spike_t = estim_spts' 
                idx_tochange = ((estim_ts>=(spike_t-pre_spike_sample*speriod))&(estim_ts<spike_t));
                sp_assoc_stimuli(idx_tochange) = estim_amps_norm(idx_tochange);
                sp_assoc_stimuli_weighted = horzcat(sp_assoc_stimuli_weighted, estim_amps_norm(idx_tochange));
            end

            ax1 = subplot(211);plot(estim_ts, estim_amps_norm,'y');hold on;
            plot(estim_ts, estim_amps_norm,'k.');hold on;
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
            sp_assoc_genSig_weighted = [];
            for spike_t = estim_spts' 
                idx_tochange = ((genSig_ts>=(spike_t-speriod))&(genSig_ts<spike_t));
                sp_assoc_genSig(idx_tochange) = genSig_vals(idx_tochange);
                sp_assoc_genSig_weighted = horzcat(sp_assoc_genSig_weighted, genSig_vals(idx_tochange)); % to be used in the next figures
            end
            %sp_assoc_genSig_weighted has the value of the genSig multiple
            %times. The count is set by the number of spikes happeing in
            %the time-window associated with that generator signal value

            ax2 = subplot(212);plot(genSig_ts, genSig_vals,'y');hold on;
            plot(genSig_ts, genSig_vals,'k.');hold on;
            plot(estim_spts, spike_genSig_vals,'r.');hold on;
            plot(genSig_ts, sp_assoc_genSig,'b*');hold on;

            xlim([0,100]);
            linkaxes([ax1,ax2],'x')

            figTitle = sprintf('%s [%s]\n trialIdx = %.2d - Extracting the spike associated stimuli/generator signal',strrep(exp_ps.exp_id,'_','.'),strrep(exp_ps.cell_id,'_','-'), trialIdx);
            suptitle(figTitle);
            saveas(gcf, [exp_ps.work_dir, fig_basename, sprintf('_F%.2d.jpeg',figIdx)]);

            %% Figure 4xx - Overlay a histogram of spike-associated generator signals.
            figIdx = 4;
            figure();

            [estim_binCounts,estim_binEdges] = histcounts(estim_amps_norm);
            estim_binCenters = (estim_binEdges(1:end-1) + estim_binEdges(2:end))/2;

            [genSig_binCounts,genSig_binEdges] = histcounts(genSig_vals);
            genSig_binCenters = (genSig_binEdges(1:end-1) + genSig_binEdges(2:end))/2;

            ax1 = subplot(321);bar(estim_binCenters, estim_binCounts,'histc');title('Normlzd Stimuli');
            ax2 = subplot(322);bar(genSig_binCenters, genSig_binCounts,'histc');title('Generator Signal');

            [sp_assoc_stimuli_binCounts,sp_assoc_stimuli_binEdges] = histcounts(sp_assoc_stimuli, estim_binEdges);
            sp_assoc_stimuli_binCenters = (sp_assoc_stimuli_binEdges(1:end-1) + sp_assoc_stimuli_binEdges(2:end))/2;
            ax3 = subplot(323);bar(sp_assoc_stimuli_binCenters, sp_assoc_stimuli_binCounts,'histc');title('Spike Associated & Normlzd Stimuli');
            
            [sp_assoc_stimuli_binCounts_weighted,~] = histcounts(sp_assoc_stimuli_weighted, estim_binEdges);
            ax5 = subplot(325);bar(sp_assoc_stimuli_binCenters, sp_assoc_stimuli_binCounts_weighted,'histc');title('Spike Associated & Weighted & Normlzd Stimuli');

            [sp_assoc_genSig_binCounts,sp_assoc_genSig_binEdges] = histcounts(sp_assoc_genSig, genSig_binEdges);
            sp_assoc_genSig_binCenters = (sp_assoc_genSig_binEdges(1:end-1) + sp_assoc_genSig_binEdges(2:end))/2;
            ax4 = subplot(324);bar(sp_assoc_genSig_binCenters, sp_assoc_genSig_binCounts,'histc');title('Spike Associated Gnerator Signal');

            [sp_assoc_genSig_binCounts_weighted,~] = histcounts(sp_assoc_genSig_weighted, genSig_binEdges);
            ax6 = subplot(326);bar(sp_assoc_genSig_binCenters, sp_assoc_genSig_binCounts_weighted,'histc');title('Spike Associated & Weighted Gnerator Signal');

            linkaxes([ax1,ax3, ax5],'x');
            linkaxes([ax2,ax4, ax6],'x');

            figTitle = sprintf('%s [%s]\n trialIdx = %.2d - Overlay a histogram of spike-associated generator signals.',strrep(exp_ps.exp_id,'_','.'),strrep(exp_ps.cell_id,'_','-'), trialIdx);
            suptitle(figTitle);
            saveas(gcf, [exp_ps.work_dir, fig_basename, sprintf('_F%.2d.jpeg',figIdx)]);

            %% Figure 5xx - [Incomplete] The probability of firing a spike vs. the generator signal
            % We would like to count the number of spikes corresponding to each value of the
            % generator signal. For this we first assign a generator signal value to
            % each spike time stamp. We then bin those spike_genSig_vals and count the
            % number of times a value falls within each bin and visualize it as a histogram.
            figIdx = 5;
            figure();

            spike_genSig_vals = interp1(genSig_ts,genSig_vals,estim_spts(estim_spts>=genSig_ts(1)));

            [FRgenSig, FRgenSig_binEdges] = histcounts(spike_genSig_vals);
            FRgenSig_binCenters = (FRgenSig_binEdges(1:end-1) + FRgenSig_binEdges(2:end))/2;

            [genSig_binCounts_tmp,~] = histcounts(genSig_vals, FRgenSig_binEdges);

            FRgenSig = FRgenSig ./ genSig_binCounts_tmp; % divide each bin by the number of generator signals there
            FRgenSig = FRgenSig .* exp_ps.stimFreq; % divide again by the sampling time to get the number of spikes per second (Hz)
            
            %bar(FRgenSig_binCenters,FRgenSig_binCounts,'histc');
            plot(FRgenSig_binCenters,FRgenSig, 'b'); hold on;
            plot(FRgenSig_binCenters,FRgenSig, 'k.');
            
            xlabel('Gen. Sig Value');
            ylabel('Spikes Rate (Hz)');

            figTitle = sprintf('%s [%s]\n trialIdx = %.2d - Number of spikes vs. the generator signal amplitude',strrep(exp_ps.exp_id,'_','.'),strrep(exp_ps.cell_id,'_','-'), trialIdx);
            suptitle(figTitle);
            saveas(gcf, [exp_ps.work_dir, fig_basename, sprintf('_F%.2d.jpeg',figIdx)]);
        end
    end
end