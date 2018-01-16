function [STA_ps, D_ps] = STA_computation(exp_ps)
    % estim_amps: electrical stimulation amplitutes for all the trials

    % This function calculates STAs, smooths the STA with cubic splining,
    % calculates parameters like location of peak & trough of STA and the
    % integration window of peak and integration window of trough
    % Develped by Nima (20180101) from basic implementation by Sudarshan(20150730)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if length(exp_ps.exp_id) > 10
        dfile = fullfile(exp_ps.data_dir,strcat(exp_ps.exp_id(1:10), exp_ps.dfile_suffix));
    else
        dfile = fullfile(exp_ps.data_dir,strcat(exp_ps.exp_id, exp_ps.dfile_suffix));        
    end
    load(dfile);
    spiketimes = eval(exp_ps.cell_id);

    stimPeriod = 1/exp_ps.stimFreq;

    if isnan(exp_ps.start_TTL) || isnan(exp_ps.stop_TTL)
        A2a = A2a(1:end, 1);% this line was not the STA_simplified
    else
        A2a = A2a(exp_ps.start_TTL:exp_ps.stop_TTL, 1);% this line was not the STA_simplified
    end
    %% ToDo: the following part of the code activated by cardinal_STA_Only_Burst should be revised
    % otherwise it expected to be highly error prone
    if exp_ps.cardinal_STA_Only_Burst
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
        Burst_End = find(diff_spiketimes > exp_ps.ISI_gap); % p.ISI_gap is what we set based on the scatter plot of pre and post ISI
        new_spike_times = [];

        for iter = 1:length(Burst_End) - 1
            if (Burst_End(iter + 1) - Burst_End(iter) > 1)
                if exp_ps.weighted_burst
                    new_spike_times_holder = repmat(spiketimes(Burst_End(iter) + 1), Burst_End(iter + 1) - Burst_End(iter), 1);
                    new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
                elseif ~exp_ps.weighted_burst
                    new_spike_times_holder = repmat(spiketimes(Burst_End(iter) + 1), 1, 1);
                    new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
                end
            end
            if (Burst_End(iter + 1) - Burst_End(iter) == 1)
                if exp_ps.singleton_spikes
                    new_spike_times_holder = repmat(spiketimes(Burst_End(iter) + 1), Burst_End(iter + 1) - Burst_End(iter), 1);
                    new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
                end
            end
        end

        if (length(spiketimes) - Burst_End(end) > 1)
            if exp_ps.weighted_burst
                new_spike_times_holder = repmat(spiketimes(Burst_End(end) + 1), length(spiketimes) - Burst_End(end), 1); % Include all the spike times after the p.last bursting event
                new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
            end
        elseif (length(spiketimes) - Burst_End(end) == 1)
            if (~exp_ps.weighted_burst || exp_ps.singleton_spikes)
                new_spike_times_holder = spiketimes(Burst_End(end) + 1);
                new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
            end
        end

        if (Burst_End(1) > 1)
            if exp_ps.weighted_burst
                new_spike_times_holder = repmat(spiketimes(1), Burst_End(1), 1); % Include all the spike times before the first bursting starts
                new_spike_times = vertcat(new_spike_times_holder, new_spike_times);
            elseif ~exp_ps.weighted_burst
                new_spike_times_holder = repmat(spiketimes(1), 1, 1); % Include all the spike times before the first bursting starts
                new_spike_times = vertcat(new_spike_times_holder, new_spike_times);
            end
        elseif (Burst_End(1) == 1)
            if (exp_ps.singleton_spikes)
                new_spike_times_holder = spiketimes(1); % Include all the spike times before the first bursting starts
                new_spike_times = vertcat(new_spike_times_holder, new_spike_times);
            end
        end
        spiketimes = new_spike_times;
    end

    %% STA computation
    tData = {};
    estim_means = [];
            
    nspikes = 0; % Number of spikes which is the spike count variable for the STAs
    nstim = 0; % Total number of stimuli vectors used for STA calculation across trials

    STA = zeros(2 * exp_ps.tKerLen, 1);
    
    if isnan(exp_ps.trials_to_use)
        trials_to_use = exp_ps.first_trial:exp_ps.last_trial;
    else
        trials_to_use = eval(exp_ps.trials_to_use);
    end
    
    flag_skip = true;
    for trialIdx = trials_to_use
        % the fist part of the following "and" is to make the code work the same as Sudas original code
        if ~isnan(exp_ps.alternate_number) % ToDo: make sure the part of the code related to flag_skip works correctly
            if (exp_ps.alternate_number > 1) && (rem(trialIdx, exp_ps.alternate_number) == 1)
                flag_skip = xor(flag_skip,true);
            else
                flag_skip = xor(flag_skip,true);
            end
        else
            flag_skip = false;
        end
        
        if ~flag_skip
            if exp_ps.NR
                estim_fname = strcat(fullfile(exp_ps.data_dir,'rexp_'), num2str(ceil(trialIdx)), '.txt');
                estim_design = importdata(estim_fname);
            else
                estim_fname = strcat(fullfile(exp_ps.data_dir,'rexp_1.txt'));
                estim_design = importdata(estim_fname);
            end

            trial_estim_amps = get_estim_amp(estim_design.textdata, exp_ps.Normalize);
            count_TTL = length(trial_estim_amps);

            trial_estim_t = A2a(((count_TTL * (trialIdx - 1)) + 1):(count_TTL * (trialIdx - 1)) + count_TTL);

            % the following line is where the shift in recorded electrical stimulation is happing
            % the source of this delay might be the equipment latency or even some cellular grounds
            trial_estim_t = trial_estim_t - exp_ps.post_wait;
            trial_estim_spt = spiketimes((spiketimes > trial_estim_t(1, 1) & spiketimes <= trial_estim_t(end, 1)), 1);

            % deletes direct RGC spikes based on lock out period
            if exp_ps.leave_out > 0
                tsp_del_loc = [];
                for esIdx = 1:length(trial_estim_t)
                    tsp_del_temp_loc = find(trial_estim_spt >= trial_estim_t(esIdx) & trial_estim_spt <= trial_estim_t(esIdx) + exp_ps.leave_out);
                    tsp_del_loc = vertcat(tsp_del_loc, tsp_del_temp_loc);
                end
                trial_estim_spt(tsp_del_loc) = [];
            end
            % Corrects spike times based on starting TTL pulse of each trial, so that all spike times will be between 0 and 100s
            trial_estim_spt = trial_estim_spt - trial_estim_t(1);
            trial_estim_t = trial_estim_t - trial_estim_t(1);

            estim_binranges = vertcat(trial_estim_t, trial_estim_t(end)+1/exp_ps.stimFreq);
            spcount_binned = histcounts(trial_estim_spt, estim_binranges);      % binned spike counts

            [trial_STA, ~, trial_nstim, trial_nspikes] = compgroupSTA_estim(trial_estim_amps, spcount_binned, exp_ps.tKerLen); % compute STA

            STA = STA + trial_STA;

            nspikes = nspikes + trial_nspikes;
            nstim = nstim + trial_nstim;

            tData(trialIdx).estim_amps = trial_estim_amps;
            tData(trialIdx).estim_ts = trial_estim_t; 
            tData(trialIdx).estim_spts = trial_estim_spt; 
            
            estim_means = [estim_means,mean(trial_estim_amps)];
            % estim_times = [estim_times; trial_estim_times];
            % the number of trial amplitutes should be 54 * 2500 = 135000
            % number of stimulus times is a subset of this total number because
            % not all stimuli are used in the STA analysis. In my opinion in the
            % outputs of the function these two should be in accordance. that
            % means the only the stimulus amplitudes used for STA should be
            % considered the same as stimulu times
        end
    end

    STA = STA' / nspikes; % Divides the Stimulus Sum/ Total Number of Spikes

    STA_t = (0.5 - exp_ps.tKerLen)* stimPeriod:stimPeriod:(exp_ps.tKerLen  - .5)* stimPeriod;
    estim_mean = mean(estim_means);

    %% splined STA
    correctedSTA = STA;
    splinedSTA_t = STA_t(1):.001:STA_t(end)+0.001;
    if exp_ps.single_pulse_activation_correction
        correctedSTA(length(STA) / 2) = estim_mean;
    end
    splinedSTA = spline(STA_t, correctedSTA, splinedSTA_t);

    %% STA significance computation
    D_ps = STA_significance(splinedSTA, estim_mean, exp_ps);
    
    %% Gather Details for later plottings
    STA_ps = struct;
    STA_ps.STA = STA;
    STA_ps.STA_t = STA_t;
    STA_ps.splinedSTA = splinedSTA;
    STA_ps.splinedSTA_t = splinedSTA_t;
    STA_ps.correctedSTA = correctedSTA;
    STA_ps.estim_mean = estim_mean;
    STA_ps.nspikes = nspikes;
    
    STA_ps.tData = tData;
    %STA_ps.estim_times = estim_times;

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

function [STA, stimulus_matrix, nstim, nspikes] = compgroupSTA_estim(estim_amp, spcount_binned, tKerLen)
    %  Inputs:
    %    estim_amp = stimulus (2D matrx: first index is time bin index, second is space)
    %    (units are the amplitude units of the stimulus, e.g. luminance or
    %    voltage)
    %    spcount_binned = column vector of spike counts in each time bin (length of 'estim_amp')
    %    tKerLen = STA kernel size. Expresed in terms of bins. Hence for a stimulation frequency of 25Hz, kernel size of 25 represents a window of 1s for which the STA is calculated.
    %    STA = Vector of stimuli preceding spikes in given bin multiplied by the number of spikes in that given bin
    %    NSP = Number of spcount_binned across all trials used for the final STA calculation
    %    NSTM = total number of stimuli vectors used for STA calculation across trials
    %    count_TTL = length of stimulus

    count_TTL = length(estim_amp);

    bins_with_spikes = find(spcount_binned); %Find which positions in the spike histogram had spikes
    spikes_before_stim = find(bins_with_spikes <= tKerLen); %Since I want to store tKerLen frames preceding each spike, I need to exclude those spikes which ocured within the
    % first tKerLen frames in a stimulus trial. So I find the location of all the spikes that occured
    % before the first tKerLen frames

    % Finds the number of spikes that happened within the last tKerLen from the end of the stimulus.

    spikes_after_stim = find(bins_with_spikes > (count_TTL - tKerLen));

    % length(location): The last spike that occured before the first tKerLen frames and find its position

    % Create a stimulus matrix of the appropriate length (by removing the spikes that occured too early or too late i.e which occured within tKerLen frames of the beginning or end of
    % stimulus block

    nstim = length(bins_with_spikes) - length(spikes_before_stim) - length(spikes_after_stim); % number of stimuli vectors to be included in STA calculation in this trial
    stimulus_matrix = zeros(nstim, tKerLen * 2); % creation of STA vector. The STA vector is tKerLen before 0s and tKerLen after 0s

    nspikes = 0; % number of spikes in each trial

    STA = zeros(2 * tKerLen, 1);

    % For loop that actually calculates the STA
    for i = 1:length(bins_with_spikes)
        if (bins_with_spikes(i) > tKerLen) && (bins_with_spikes(i) < (count_TTL - tKerLen)) % If the spike occurs after the first tKerLen frames in a trial or before the last tKerLen frames in a trial include the stimuli for the STA calculation
            stimulus_matrix(i, :) = estim_amp(:, (bins_with_spikes(i)) + (1 - tKerLen):(bins_with_spikes(i)) + tKerLen)'; % Stores the stimuli used for the STA calculation
            STA = STA + (stimulus_matrix(i, :) .* spcount_binned(bins_with_spikes(i)))';
            % The STA vector adds up all the stimuli occuring before a spike. The stimuli are multiplied by the number of spikes they cause before they are used for the STA calculation
            nspikes = nspikes + spcount_binned(bins_with_spikes(i)); % Update the number of spikes used for the STA calculation        
        end
    end
end