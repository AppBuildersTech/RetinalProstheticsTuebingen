function [STA_ps, exp_ps] = STA_computation(exp_ps)
    % This function calculates STAs, smooths the STA with cubic splining,
    % calculates parameters like location of peak & trough of STA and the
    % integration window of peak and integration window of trough
    % Develped by Nima (20180101) from basic implementation by Sudarshan(20150730)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Get Data
    if length(exp_ps.exp_id) > 10
        dfile = fullfile(exp_ps.data_dir,strcat(exp_ps.exp_id(1:10), exp_ps.dfile_suffix));
    else
        dfile = fullfile(exp_ps.data_dir,strcat(exp_ps.exp_id, exp_ps.dfile_suffix));        
    end
    stimPeriod = 1/exp_ps.stimFreq;
    
    trial_nsample = exp_ps.trial_length_in_secs *exp_ps.stimFreq;

    load(dfile);
    spiketimes = eval(exp_ps.cell_id);

    estim_t = A2a((1+exp_ps.cut_TTL_head):(end-exp_ps.cut_TTL_tail), 1); % stimulus timings
    
    experiment_nsample = length(estim_t);

    % Get Stimulus amplitudes
    if ~isnan(exp_ps.mat_stim)
        estim_fname = strcat(fullfile(exp_ps.data_dir,exp_ps.mat_stim), '.mat');
        tmp_estim_amps = importdata(estim_fname);
        estim_amps = zeros(experiment_nsample,1);
        for trialIdx = 1:size(tmp_estim_amps,1)
            estim_amps(((trialIdx-1)*trial_nsample+1):trialIdx*trial_nsample,1) = tmp_estim_amps(trialIdx, :);
        end
    else
        estim_amps = zeros(experiment_nsample,1);
        for trialIdx = 1:experiment_nsample/trial_nsample
            if ~exp_ps.single_stim
                estim_fname = strcat(fullfile(exp_ps.data_dir,'rexp_'), num2str(ceil(trialIdx)), '.txt');
                temp_estim = importdata(estim_fname);
            else
                estim_fname = strcat(fullfile(exp_ps.data_dir,'rexp_1.txt'));
                temp_estim = importdata(estim_fname);
            end
            estim_amps(((trialIdx-1)*trial_nsample+1):trialIdx*trial_nsample,1) = get_estim_amp(temp_estim.textdata, exp_ps.Normalize);
        end
    end
    
    %% Exclude spikes in between trials
    for tIdx = 2:length(estim_t)
        if (estim_t(tIdx) - estim_t(tIdx-1))>stimPeriod
            spiketimes(spiketimes>estim_t(tIdx-1) & spiketimes<estim_t(tIdx)) = [];
        end
    end
    
    %% Exclude direct RGC spikes based on lock out period
    if exp_ps.leave_out > 0
        tsp_del_loc = [];
        for esIdx = 1:length(estim_t)
            tsp_del_temp_loc = find(spiketimes >= estim_t(esIdx) & spiketimes <= estim_t(esIdx) + exp_ps.leave_out);
            tsp_del_loc = vertcat(tsp_del_loc, tsp_del_temp_loc);
        end
        spiketimes(tsp_del_loc) = [];
    end
        
    estim_t = estim_t - exp_ps.post_wait;
   
    %% Exclude spikes that occure before exp_ps.tKerLen
    spiketimes(spiketimes<exp_ps.tKerLen*stimPeriod) = [];
    
    %% ToDo: the following part of the code activated by BTA should be revised
    % otherwise it expected to be highly error prone
    if exp_ps.BTA
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
    spcount = horzcat(0,histcounts(spiketimes,estim_t))';
    spcount = spcount(exp_ps.tKerLen:end); % a critical step        
    nsp = sum(spcount); % Number of spikes which is the spike count variable for the STAs

    % the following line is where the shift in recorded electrical stimulation is happing
    % the source of this delay might be the equipment latency or even some cellular grounds

    raw_stim_ensem = getStimSegments(estim_amps,exp_ps.tKerLen); % Convert stimulus to a matrix where each row is one segment of stim
    
    stim_len = length(estim_t) - exp_ps.tKerLen; 
    % Compute raw mean and covariance
    RawMu = mean(raw_stim_ensem)';
    RawCov = (raw_stim_ensem'*raw_stim_ensem)/(stim_len-1)-RawMu*RawMu'*stim_len/(stim_len-1);

    % Compute spike-triggered mean and covariance
    spikeIds = find(spcount>0);
    spikeIds((spikeIds+exp_ps.tKerLen)>(experiment_nsample-exp_ps.tKerLen+1)) = [];
    
    spike_assoc_ensemble = raw_stim_ensem(spikeIds,:);
    spike_assoc_ensemble_spcount = spcount(spikeIds);
    
    STA = (spike_assoc_ensemble_spcount'*spike_assoc_ensemble)'/nsp;
    STA_t = ((0.5 - exp_ps.tKerLen)* stimPeriod:stimPeriod:0)';% we want to show bin centers

    STC = spike_assoc_ensemble'*(spike_assoc_ensemble.*repmat(spike_assoc_ensemble_spcount,1,exp_ps.tKerLen))/(nsp-1) - STA*STA'*nsp/(nsp-1);
    
    postSTA = (spike_assoc_ensemble_spcount'*raw_stim_ensem(spikeIds+exp_ps.tKerLen,:))'/nsp;%will be used for significance check and visualization purposes
    postSTA_t = (STA_t(end)+stimPeriod:stimPeriod:(exp_ps.tKerLen  - .5)* stimPeriod)';
    
    estim_mean = mean(estim_amps);
    estim_std = std(estim_amps);% How much variance on average we have

    %% splined STA and STA significance computation
    D_ps = STA_significance(STA, postSTA, estim_mean, exp_ps);
    
    % when the exp_ps.single_pulse_activation_correction is set then the
    % returned STA will be zero time point correceted. If not then the program
    % will still correct the zero time point (to clear the large negative
    % deflection) and calculate D points but will return the non-corrected
    % STA as for the suplined STA. This point should be made known to the user
    if exp_ps.single_pulse_activation_correction % NOTE THIS IF and the above comments
        STA(end) = estim_mean;
    end
    
    %% Gather Details for later plottings
    STA_ps = struct;
    
    STA_ps.STA = STA;
    STA_ps.STA_t = STA_t;
    
    STA_ps.postSTA = postSTA;
    STA_ps.postSTA_t = postSTA_t;
    
    STA_ps.D_ps = D_ps;
    
    exp_ps.estim_mean = estim_mean;
    exp_ps.estim_std = estim_std;
    exp_ps.nspikes = nsp;
    
    exp_ps.tData = struct;
    exp_ps.tData.estim_amps = estim_amps;
    exp_ps.tData.estim_ts = estim_t; 
    exp_ps.tData.estim_spts = spiketimes;
    
    exp_ps.tData.raw_stim_ensem = raw_stim_ensem;
    exp_ps.tData.spike_assoc_ensemble = spike_assoc_ensemble;
    exp_ps.tData.spike_assoc_ensemble_spcount = spike_assoc_ensemble_spcount;

end

function raw_stim_ensem = getStimSegments(estim,Kw)
    %% return segmented stimuli, where each segment is Kw length and only 1 sample apart from the next segment
    T = size(estim,1);
    raw_stim_ensem = zeros(T-Kw+1,Kw);
    for xIdx = 1:T-Kw+1
        raw_stim_ensem(xIdx,:) = estim(xIdx:(xIdx+Kw-1),1);
    end
end

function trial_estim_amps = get_estim_amp(textdata, Normalize)
    trial_estim_amps = [];
    % formats stimuli from text file into a vector
    for i = 8:2:length(textdata)%p.last % 8 is the starting line in the estim_amps(trialIdx,:) textfile
        x = str2double(textdata{i, 1});
        trial_estim_amps = horzcat(trial_estim_amps, x);
    end
    % Normalises the trial_estim_amps if specified by user
    if Normalize
        trial_estim_amps = (trial_estim_amps - mean(trial_estim_amps)) / std(trial_estim_amps); %stim = Stimulus;% Sudsa Modified
    end
end