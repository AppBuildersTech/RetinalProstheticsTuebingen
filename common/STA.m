function [estim_amps_all stimulus_matrix V_excite V_inhibit Cov_Matrix Cov_Matrix_raw_stimulus] = STA(p)

% estim_amps_all: electrical stimulation amplitutes for all the trials

% This function calculates STAs, smooths the STA with cubic splining,
% calculates parameters like location of peak & trough of STA and the
% integration window of peak and integration window of trough
% In this function I call 1) compgroupSTA_estim 2) STA_parameters_significance 3) STC_excitatory_parameters_significance 4) STC_inhibitory_parameters_significance
% Created by Sudsa (20150730)
% Improved by Nima (20180101)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(p.data_dir,strcat(p.exp_id, p.dfile_suffix)));
spiketimes = eval(p.cell_id);

if isnan(p.start_TTL) || isnan(p.stop_TTL)
    A2a = A2a(1:end, 1);% this line was not the STA_simplified
else
    A2a = A2a(p.start_TTL:p.stop_TTL, 1);% this line was not the STA_simplified
end
%% ToDo: the following part of the code activated by cardinal_STA_Only_Burst should be revised
% otherwise it expected to be highly error prone
if p.cardinal_STA_Only_Burst
    % fsta spike times are recalculated so that the original spike train is burst corrected.
    % The corrected spike times can be weighted.
    % Singleton spikes can be included.
    % spiketimes(Burst_End(iter)+1) gives the first spike right after a burst.
    % Then the subsequesnt spiketimes in a burst are replaced by the first spike of the burst.
    % The number of spikes to be replaced is given by the line
    % Burst_End(iter+1)-Burst_End(iter). This is basically the length of the burst.     
    % The entire operation is being done in a for loop.
    % So in each iteration the spikes are concatenated to get the updated spiketimes
        
    diff_spiketimes = diff(spiketimes); % diff_spiketimes
    Burst_End = find(diff_spiketimes > p.ISI_gap); % p.ISI_gap is what we set based on the scatter plot of pre and post ISI
    new_spike_times = [];
 
    for iter = 1:length(Burst_End) - 1
        if (Burst_End(iter + 1) - Burst_End(iter) > 1)
            if p.weighted_burst
                new_spike_times_holder = repmat(spiketimes(Burst_End(iter) + 1), Burst_End(iter + 1) - Burst_End(iter), 1);
                new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
            elseif ~p.weighted_burst
                new_spike_times_holder = repmat(spiketimes(Burst_End(iter) + 1), 1, 1);
                new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
            end
        end
        if (Burst_End(iter + 1) - Burst_End(iter) == 1)
            if p.singleton_spikes
                new_spike_times_holder = repmat(spiketimes(Burst_End(iter) + 1), Burst_End(iter + 1) - Burst_End(iter), 1);
                new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
            end
        end
    end

    if (length(spiketimes) - Burst_End(end) > 1)
        if p.weighted_burst
            new_spike_times_holder = repmat(spiketimes(Burst_End(end) + 1), length(spiketimes) - Burst_End(end), 1); % Include all the spike times after the p.last bursting event
            new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
        end
    elseif (length(spiketimes) - Burst_End(end) == 1)
        if (~p.weighted_burst || p.singleton_spikes)
            new_spike_times_holder = spiketimes(Burst_End(end) + 1);
            new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
        end
    end
 
    if (Burst_End(1) > 1)
        if p.weighted_burst
            new_spike_times_holder = repmat(spiketimes(1), Burst_End(1), 1); % Include all the spike times before the first bursting starts
            new_spike_times = vertcat(new_spike_times_holder, new_spike_times);
        elseif ~p.weighted_burst
            new_spike_times_holder = repmat(spiketimes(1), 1, 1); % Include all the spike times before the first bursting starts
            new_spike_times = vertcat(new_spike_times_holder, new_spike_times);
        end
    elseif (Burst_End(1) == 1)
        if (p.singleton_spikes)
            new_spike_times_holder = spiketimes(1); % Include all the spike times before the first bursting starts
            new_spike_times = vertcat(new_spike_times_holder, new_spike_times);
        end
    end
    spiketimes = new_spike_times;
end

%% calculating spiketimes within each trial for a cell
estim_amps_all = []; % electrical stimulation amplitutes for all the trials
stimulus_matrix = [];

NSP = 0; % Number of spikes which is the spike count variable for the STAs
NSTM = 0; % Total number of stimuli vectors used for STA calculation across trials


p.skip_cycle = p.alternate_number * 2; % this line was not the STA_simplified

STA = zeros(2 * p.tKerLen, 1);
flag_skip = true;

tmp_NSP = 0;
tmp_inner = 0;

if isnan(p.trials_to_use)
    trials_to_use = p.first_trial:p.last_tiral;
else
    trials_to_use = p.trials_to_use;
end
    
% Goes through trials
for trialIdx = trials_to_use
    % the fist part of the following "and" is to make the code work the same as Sudas original code
    if isnan(p.trials_to_use) && p.alternate
        if (p.alternate_number > 1) && (rem(trialIdx, p.alternate_number) == 1)
            flag_skip = xor(flag_skip,true);
        else
            flag_skip = xor(flag_skip,true);
        end
    else
        flag_skip = false;
    end
    if ~flag_skip
        if p.NR
            estim_fname = strcat(fullfile(p.data_dir,'rexp_'), num2str(ceil(trialIdx)), '.txt');
            estim_design = importdata(estim_fname);
        end
        estim_amp = get_estim_amp(estim_design.textdata, p.Normalize);
        count_TTL = length(estim_amp);
        
        % the fist part of the following "or" is to make the code work the same as Sudas original code
        % note that estim_timings will hold the start time of the applied stimulus
        if isnan(p.trials_to_use) || (p.leaveit == 1)
            estim_timings = A2a(((count_TTL * (trialIdx - 1)) + 1):(count_TTL * (trialIdx - 1)) + count_TTL);
        elseif ~isnan(p.trials_to_use) && (p.leaveit == 2)
            if rem(trialIdx, p.alternate_number) == 0
                estim_timings = A2a((floor(trialIdx / p.skip_cycle) * (p.alternate_number * (count_TTL + p.skip))) + (count_TTL * (p.alternate_number - 1)) + 1:(floor(trialIdx / p.skip_cycle) * (p.alternate_number * ...
                (count_TTL + p.skip))) + (count_TTL * p.alternate_number));
            else
                estim_timings = A2a((floor(trialIdx / p.skip_cycle) * (p.alternate_number * (count_TTL + p.skip))) + (count_TTL * (rem(trialIdx, p.alternate_number) - 1)) + 1:(floor(trialIdx / p.skip_cycle) * (p.alternate_number * ...
                (count_TTL + p.skip))) + (count_TTL * rem(trialIdx, p.alternate_number)));
            end
        end
        % the following line is where the shift in recorded electrical stimulation is happing
        % the source of this delay might be the equipment latency or even some cellular grounds
        estim_timings = estim_timings - p.post_wait;
        estim_spiketimes = spiketimes((spiketimes > estim_timings(1, 1) & spiketimes <= estim_timings(end, 1)), 1);

        % deletes direct RGC spikes based on lock out period
        if p.leave_out > 0
            tsp_del_loc = [];
            for esIdx = 1:length(estim_timings)
                tsp_del_temp_loc = find(estim_spiketimes >= estim_timings(esIdx) & estim_spiketimes <= estim_timings(esIdx) + p.leave_out);
                tsp_del_loc = vertcat(tsp_del_loc, tsp_del_temp_loc);
            end
            estim_spiketimes(tsp_del_loc) = [];
        end
        % Corrects spike times based on starting TTL pulse of each trial, so that all spike times will be between 0 and 100s
        estim_spiketimes = estim_spiketimes - estim_timings(1);
        estim_timings = estim_timings - estim_timings(1);

        estim_binranges = vertcat(estim_timings, estim_timings(end)+1/p.stimFreq);
        spcount_binned = histcounts(estim_spiketimes, estim_binranges);      % binned spike counts

        [STA, stimulus_matrix, nstim , nspikes, Stm_matrix] = compgroupSTA_estim(estim_amp, spcount_binned, p.tKerLen, STA); % compute STA
        
        NSP = NSP + nspikes;
        NSTM = NSTM + nstim;
        
        estim_amps_all = [estim_amps_all estim_amp]; 
        stimulus_matrix = [stimulus_matrix; Stm_matrix];
        
        tmp_NSP = tmp_NSP + length(estim_spiketimes);
        tmp_inner = tmp_inner + length(estim_timings);
    end
end

%% Plots begin from here
line_thickness = 2;

if p.Normalize == 1
    plt_ylim = [- 1, 1];
else
	plt_ylim = [- 1300, -300];
end

STA = STA / NSP; % Divides the Stimulus Sum/ Total Number of Spikes
y = 1 * STA;
frame = 1 / (count_TTL / p.trial_length_in_secs);
x = (- (p.tKerLen)) * (frame) + .5 * frame:frame:p.tKerLen * frame - .5 * frame;
h = plot(x, y);
set(h, 'LineWidth', line_thickness)
baseline = (1 * mean(estim_amp)) * ones(2 * p.tKerLen, 1);
hold on
plot(x, baseline, 'k')
set(gcf, 'color', 'w');
zeromarker = zeros(length(plt_ylim(1):100:plt_ylim(2)));
hold on
plot(zeromarker, plt_ylim(1):100:plt_ylim(2), 'k');
ylim([plt_ylim(1) plt_ylim(2)])
hold on
legend(['STA:spks = ', int2str(NSP)], ['Mean Stimulus : ', int2str(p.stimFreq), 'Hz stm. frq. : 35% var '])
%%
fig_names = strcat(p.names, p.cell_id);
total_time = p.trial_length_in_secs * length(trials_to_use);
Frequency = NSP / (length(trials_to_use) * p.trial_length_in_secs);
xlabel({'Time (sec)'})
ylabel('Mean Stimulus (mV)')
mousetype = ' C57Bl/6 (wt) ';
figure_heading = strcat(p.year, p.cell_id, mousetype, ' FR = ', num2str(Frequency), ' Hz: ', num2str(total_time), ' sec : epiretinal', ' lock out is ', num2str(p.leave_out), 's', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes));
title(figure_heading)
set(gcf, 'name', fig_names)
fname = strcat(p.exp_id, int2str(p.stimFreq), 'Hz_', num2str(p.cont), '_', p.mean_V, '\', p.cell_id)
fname = strcat(p.work_dir, fname, '\');
if ~ exist(fname, 'dir'), mkdir(fname); end
 
    ext = '.fig';
 
    if (p.Normalize == 0)
     
        saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'not_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), ' WB is ', num2str(p.weighted_burst), ' Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
        set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
     
        saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'not_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), ' WB is ', num2str(p.weighted_burst), ' Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
     
        % close all
        figure
        %%
        y = 1 * STA;
        frame = 1 / (count_TTL / p.trial_length_in_secs);
     
        x = (- (p.tKerLen)) * (frame) + .5 * frame:frame:p.tKerLen * frame - .5 * frame;
        h = plot(x, y);
        set(gcf, 'color', 'w');
     
        saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'not_norm_for_NL', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
        set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
     
        saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'not_norm_for_NL', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
     
        %% error bars
     
        figure
        y = 1 * STA;
        STA_error = std(y(length(y) / 2 + 1:end)) / sqrt(p.tKerLen);
        STA_error = STA_error * ones(length(y), 1);
        dim = max(size(stimulus_matrix));
        for ijk = 1:2 * p.tKerLen
            errorbar_sta(ijk) = std(stimulus_matrix(:, ijk)) / sqrt(dim);
        end
     
        frame = 1 / (count_TTL / p.trial_length_in_secs);
        x = (- (p.tKerLen)) * (frame) + .5 * frame:frame:p.tKerLen * frame - .5 * frame;
        h = errorbar(x, y, errorbar_sta);
        set(h, 'LineWidth', line_thickness)
        baseline = (1 * mean(estim_amp)) * ones(2 * p.tKerLen, 1);
        hold on
     
        errorbar(x, baseline, STA_error, 'k')
        set(gcf, 'color', 'w');
        zeromarker = zeros(length(plt_ylim(1):100:plt_ylim(2)));
        hold on
        plot(zeromarker, plt_ylim(1):100:plt_ylim(2), 'k');
        ylim([plt_ylim(1) plt_ylim(2)])
        hold on
     
        legend(['STA:spks = ', int2str(NSP)], ['Mean Stimulus : ', int2str(p.stimFreq), 'Hz stm. frq. : 35% var '])
        total_time = p.trial_length_in_secs * length(trials_to_use);
        Frequency = NSP / (length(trials_to_use) * p.trial_length_in_secs);
        xlabel({'Time (sec)'})
        ylabel('Mean Stimulus (mV)')
        mousetype = ' C57Bl/6 (wt) ';
     
        figure_heading = strcat(p.year, p.cell_id, mousetype, ' FR = ', num2str(Frequency), ' Hz: ', num2str(total_time), ' sec : epiretinal', ' lock out is ', num2str(p.leave_out), 's', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes));
        title(figure_heading)
        set(gcf, 'name', fig_names)
        if ~ exist(fname, 'dir'), mkdir(fname); end
            ext = '.fig';
            saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'error_bars_not_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
            set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
            saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'error_bars_not_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
         
        else
         
            set(h, 'LineWidth', line_thickness)
            baseline = (1 * mean(estim_amp)) * ones(2 * p.tKerLen, 1);
            hold on
            plot(x, baseline, 'k')
            set(gcf, 'color', 'w');
            zeromarker = zeros(length(plt_ylim(1):.05:plt_ylim(2)));
            hold on
            plot(zeromarker, plt_ylim(1):.05:plt_ylim(2), 'k');
            ylim([plt_ylim(1) plt_ylim(2)])
            hold on
         
            legend(['STA:spks = ', int2str(NSP)], ['Mean Stimulus : ', int2str(p.stimFreq), 'Hz stm. frq. : 35% var '])
         
            %%
            total_time = p.trial_length_in_secs * length(trials_to_use);
            Frequency = NSP / (length(p.first_trial:p.last_tiral) * p.trial_length_in_secs);
            xlabel({'Time (sec)'})
            ylabel('Mean Stimulus (mV)')
            mousetype = ' C57Bl/6 (wt) ';
            set(gcf, 'color', 'w');
            saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'p.Normalised_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
            set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
         
            saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'p.Normalised_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
            figure
            y = 1 * STA;
            frame = 1 / (count_TTL / p.trial_length_in_secs);
            x = (- (p.tKerLen)) * (frame) + .5 * frame:frame:p.tKerLen * frame - .5 * frame;
            h = plot(x, y);
            set(gcf, 'color', 'w');
            saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'p.Normalised_for_NL', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
            set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
            saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'p.Normalised_for_NL', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
         
            %% error bars
            y = 1 * STA;
            STA_error = std(y(length(y) / 2 + 1:end)) / sqrt(p.tKerLen);
            STA_error = STA_error * ones(length(y), 1);
            frame = 1 / (count_TTL / p.trial_length_in_secs);
            x = (- (p.tKerLen)) * (frame) + .5 * frame:frame:p.tKerLen * frame - .5 * frame;
            h = plot(x, y);
            set(h, 'LineWidth', line_thickness)
            baseline = (1 * mean(estim_amp)) * ones(2 * p.tKerLen, 1);
            hold on
            plot(zeromarker, plt_ylim(1):.05:plt_ylim(2), 'k');
            errorbar(x, baseline, STA_error, 'k')
            zeromarker = zeros(length(plt_ylim(1):100:plt_ylim(2)));
            hold on
            set(gcf, 'color', 'w');
            ylim([plt_ylim(1) plt_ylim(2)])
            hold on
            legend(['STA:spks = ', int2str(NSP)], ['Mean Stimulus : ', int2str(p.stimFreq), 'Hz stm. frq. : 35% var '])
         
            %%
            total_time = p.trial_length_in_secs * length(trials_to_use);
            Frequency = NSP / (length(p.first_trial:p.last_tiral) * p.trial_length_in_secs);
            xlabel({'Time (sec)'})
            ylabel('Mean Stimulus (mV)')
            mousetype = ' C57Bl/6 (wt) ';
         
            figure_heading = strcat(p.year, p.cell_id, mousetype, ' FR = ', num2str(Frequency), ' Hz: ', num2str(total_time), ' sec : epiretinal', ' lock out is ', num2str(p.leave_out), 's', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes));
            title(figure_heading)
            set(gcf, 'name', fig_names)
            if ~ exist(fname, 'dir'), mkdir(fname); end
                ext = '.fig';
                saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'error_bars_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
                set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
                saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'error_bars_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
            end
            % STA parameter calculation
            STA_parameters_significance(Frequency, fname, frame, STA, estim_amp, p);
            % STC
            if p.STC_Analysis
                [V_inhibit V_excite Cov_Matrix Cov_Matrix_raw_stimulus mu sd] = STC_excitatory_parameters_significance(p.tKerLen, frame, STA, estim_amps_all, stimulus_matrix, p.cell_id, p.year, p.cardinal_STA_Only_Burst);
                %  [V_inhibit Cov_Matrix Cov_Matrix_raw_stimulus estim_amps_all] = STC_inhibitory_parameters_significance(p.tKerLen, frame, STA, estim_amps_all, stimulus_matrix, p.cell_id, p.year, p.cardinal_STA_Only_Burst, mu, sd);
            else
                V_excite = 0;
                V_inhibit = 0;
                Cov_Matrix = 0;
                Cov_Matrix_raw_stimulus = 0;
            end
end
        
function estim_amp = get_estim_amp(textdata, Normalize)
    estim_amp = [];
    % formats stimuli from text file into a vector
    for i = 8:2:length(textdata)%p.last % 8 is the starting line in the estim_amp textfile
        x = str2double(textdata{i, 1});
        estim_amp = horzcat(estim_amp, x);
    end
    % Normalises the estim_amp if specified by user
    if Normalize
        estim_amp = (estim_amp - mean(estim_amp)) / std(estim_amp); %stim = Stimulus;% Sudsa Modified
    end
end