function [estim_amps estim_times V_excite V_inhibit Cov_Matrix Cov_Matrix_raw_stimulus] = STA(p)

% estim_amps: electrical stimulation amplitutes for all the trials

% This function calculates STAs, smooths the STA with cubic splining,
% calculates parameters like location of peak & trough of STA and the
% integration window of peak and integration window of trough
% In this function I call 1) compgroupSTA_estim 2) STA_parameters_significance 3) STC_excitatory_parameters_significance 4) STC_inhibitory_parameters_significance
% Created by Sudsa (20150730)
% Improved by Nima (20180101)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(p.data_dir,strcat(p.exp_id, p.dfile_suffix)));
spiketimes = eval(p.cell_id);

stimPeriod = 1/p.stimFreq;

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

%% STA computation
estim_amps = []; % electrical stimulation amplitutes for all the trials
estim_times = [];

nspikes = 0; % Number of spikes which is the spike count variable for the STAs
nstim = 0; % Total number of stimuli vectors used for STA calculation across trials

p.skip_cycle = p.alternate_number * 2; % this line was not the STA_simplified

STA = zeros(2 * p.tKerLen, 1);

flag_skip = true;

if isnan(p.trials_to_use)
    trials_to_use = p.first_trial:p.last_trial;
else
    trials_to_use = p.trials_to_use;
end
    
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
        trial_estim_amps = get_estim_amp(estim_design.textdata, p.Normalize);
        count_TTL = length(trial_estim_amps);
        
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

        [trial_STA, trial_estim_times, trial_nstim, trial_nspikes] = compgroupSTA_estim(trial_estim_amps, spcount_binned, p.tKerLen); % compute STA
        
        STA = STA + trial_STA;
        
        nspikes = nspikes + trial_nspikes;
        nstim = nstim + trial_nstim;
        
        estim_amps = [estim_amps; trial_estim_amps]; 
        estim_times = [estim_times; trial_estim_times];
        % the number of trial amplitutes should be 54 * 2500 = 135000
        % number of stimulus times is a subset of this total number because
        % not all stimuli are used in the STA analysis. In my opinion in the
        % outputs of the function these two should be in accordance. that
        % means the only the stimulus amplitudes used for STA should be
        % considered the same as stimulu times
    end
end

STA = STA / nspikes; % Divides the Stimulus Sum/ Total Number of Spikes

%% Plots begin from here
line_thickness = 2;

STA_t = (0.5 - p.tKerLen)* stimPeriod:stimPeriod:(p.tKerLen  - .5)* stimPeriod;
estim_mean = ceil(mean(mean(estim_amps)));
estim_meanline = estim_mean * ones(2 * p.tKerLen, 1);
total_trial_time = p.trial_length_in_secs * length(trials_to_use);
mean_FR = nspikes / total_trial_time;

if p.Normalize == 1
    plt_ylim = [-1, 1];
else
	plt_ylim = [-1300, -300];
end

fig_basename = sprintf('[%s]_%s_%dto%d_FR=%2.3fHz_leave_out=%0.2fs_cSOB=%d_WB=%d_SS=%d',...
    p.cell_id,p.year,p.first_trial,p.last_trial,mean_FR,p.leave_out,p.cardinal_STA_Only_Burst,p.weighted_burst,p.singleton_spikes);

%% Plot 1: STA
figure
set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
set(gcf, 'color', 'w');

plot(STA_t, STA,'LineWidth', line_thickness);hold on;
plot(STA_t, estim_meanline, 'k');

yaxis_line = zeros(length(plt_ylim(1):100:plt_ylim(2)));
plot(yaxis_line, plt_ylim(1):100:plt_ylim(2), 'k');

title(fig_basename, 'Interpreter', 'none')
ylim([plt_ylim(1) plt_ylim(2)])

saveas(gcf, [p.work_dir, fig_basename,'_STA plot.fig']);
saveas(gcf, [p.work_dir, fig_basename, '_STA plot.jpeg']);

%% Plot 2: STA with error bars - ToDo: write STA with errorbars again IMHO code was either extremely complex written or not correct for their intented purpose
% STA_error = std(STA(length(STA) / 2 + 1:end)) / sqrt(p.tKerLen);
% STA_error = STA_error * ones(length(STA), 1);
% errorbar(STA_t, estim_meanline, STA_error, 'k')

%% STA significance computation
STA_parameters_significance(STA, STA_t, trial_estim_amps, p);
%STA_parameters_significance(single_pulse_corrected_STA, csSTA, STA_t, csSTA_t, trial_estim_amps, p);
% this the original implementation. notice how the last trials stimulus amplitude is passed to this function! 

% STC
if p.STC_Analysis
    [V_inhibit V_excite Cov_Matrix Cov_Matrix_raw_stimulus mu sd] = STC_excitatory_parameters_significance(p.tKerLen, stimPeriod, STA, estim_amps, estim_times, p.cell_id, p.year, p.cardinal_STA_Only_Burst);
    %  [V_inhibit Cov_Matrix Cov_Matrix_raw_stimulus estim_amps] = STC_inhibitory_parameters_significance(p.tKerLen, stimPeriod, STA, estim_amps, estim_times, p.cell_id, p.year, p.cardinal_STA_Only_Burst, mu, sd);
else
    V_excite = 0;
    V_inhibit = 0;
    Cov_Matrix = 0;
    Cov_Matrix_raw_stimulus = 0;
end

end
        
function trial_estim_amps = get_estim_amp(textdata, Normalize)
    trial_estim_amps = [];
    % formats stimuli from text file into a vector
    for i = 8:2:length(textdata)%p.last % 8 is the starting line in the trial_estim_amps textfile
        x = str2double(textdata{i, 1});
        trial_estim_amps = horzcat(trial_estim_amps, x);
    end
    % Normalises the trial_estim_amps if specified by user
    if Normalize
        trial_estim_amps = (trial_estim_amps - mean(trial_estim_amps)) / std(trial_estim_amps); %stim = Stimulus;% Sudsa Modified
    end
end