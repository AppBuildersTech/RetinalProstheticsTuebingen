function [STA, stimulus_matrix, nstim, nspikes] = compgroupSTA_estim(estim_amp, spcount_binned, tKerLen)
% COMPutes GROUP (multiple repeats of the stimulus) Spike Triggered Average for Electrical STIMulation data.
%
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

% This code is mirrored compgroupSTA.m for visual stimuli.

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
