function STA_parameters_significance(STA, STA_t, correctedSTA, splinedSTA, splinedSTA_t, trial_estim_amps, estim_meanline, p)

% This function calculates calculates parameters like location of peak & trough of STA and the
% integration window of peak and integration window of trough
% fig4 is the cubic splined STA. Marked with black circles, are the width of
% the integration window of peak and width of integration windown of trough
% STA calls this function

% estim_amps is the stimulus of the last stimulus trial. It is used calculate estim_meanline mean and variance which later are used for significance calculation later.

stimPeriod = 1/p.stimFreq;

if p.Normalize == 1
    plt_ylim = [-1, 1];
else
    plt_ylim = [-1300, -300];
end

%% Significance of the STA peak and trough

STA_mean = mean(STA(length(STA) / 2 + 1:end));
STA_std = std(STA(length(STA) / 2 + 1:end));

alpha_value = (1 - nthroot(.95, p.tKerLen))/2;

[~, STA_peak_idx] = max(abs(STA(1:length(STA) / 2)));
[~, STA_trough_idx] = min(abs(STA(1:length(STA) / 2)));
STA_peak_val = STA(STA_peak_idx);
STA_trough_val = STA(STA_trough_idx);

[STA_peak_issig, ~] = ztest(STA_peak_val, STA_mean, STA_std, 'Alpha',alpha_value); % h = 1 means an outlier, and we want an outlier
[STA_through_issig, ~] = ztest(STA_trough_val, STA_mean, STA_std, 'Alpha', alpha_value);

%% Location of the peak and the trough and their values

% [splinedSTA_peak_val, splinedSTA_peak_idx] = min(splinedSTA(1:(length(splinedSTA) / 2)));
% [splinedSTA_trough_val, splinedSTA_trough_idx] = max(splinedSTA(1:length(splinedSTA) / 2));

[~, splinedSTA_peak_idx] = max(abs(splinedSTA(1:length(splinedSTA) / 2)));
[~, splinedSTA_trough_idx] = min(abs(splinedSTA(1:length(splinedSTA) / 2)));
splinedSTA_peak_val = splinedSTA(splinedSTA_peak_idx);
splinedSTA_trough_val = splinedSTA(splinedSTA_trough_idx);

if ~ p.vstim
    trough_location_time = abs(splinedSTA_t(1) + (splinedSTA_trough_idx -1)*.001);
    peak_location_time = abs(splinedSTA_t(1) + (splinedSTA_peak_idx  - 1)*.001);
elseif p.vstim
    trough_location_time = abs(splinedSTA_t(1) + (splinedSTA_peak_idx -1)*.001);
    peak_location_time = abs(splinedSTA_t(1) + (splinedSTA_trough_idx  - 1)*.001);
end

%% integration time calculation peak

if ~STA_peak_issig; splinedSTA_peak_idx = - 1000; end % if peak is not an outlier
if ~STA_through_issig; splinedSTA_trough_idx = - 1000; end

if (splinedSTA_peak_idx < splinedSTA_trough_idx)
        
    for i_zerocrossing = splinedSTA_trough_idx:-1:1%(splinedSTA(i_zerocrossing) > mean(trial_estim_amps))
        if ~(splinedSTA(i_zerocrossing) > mean(trial_estim_amps)); break; end
    end

    for i_sig = splinedSTA_trough_idx:-1:1
        if ~((splinedSTA(i_sig) > mean(trial_estim_amps))&& ztest(splinedSTA(i_sig), STA_mean, STA_std, 'Alpha',alpha_value)); break; end
    end

    plot(STA_t(1) - .001 + (i_zerocrossing * .001), splinedSTA(i_zerocrossing), 'ok', 'LineWidth', 2)
    hold on
    plot(STA_t(1) - .001 + (i_sig * .001), splinedSTA(i_sig), 'or', 'LineWidth', 2)
    
    for j_zerocrossing = splinedSTA_trough_idx:(length(splinedSTA) / 2)+1 % BUG: the previous code counted always one more than the 
        if ~(splinedSTA(j_zerocrossing) >= mean(trial_estim_amps)); break; end
    end
    
    for j_sig = splinedSTA_trough_idx:(length(splinedSTA) / 2)+1 % BUG: the previous code counted always one more than the 
        if ~((splinedSTA(j_sig) > mean(trial_estim_amps))&& ztest(splinedSTA(j_sig), STA_mean, STA_std, 'Alpha',alpha_value)); break; end
    end
    
    plot(STA_t(1) - .001 + (j_zerocrossing * .001), splinedSTA(j_zerocrossing), 'ok', 'LineWidth', 2)
    hold on
    plot(STA_t(1) - .001 + (j_sig * .001), splinedSTA(j_sig), 'or', 'LineWidth', 2)
    
    if ~ p.vstim
        trough_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
        trough_Integration_time_sig = (j_sig - i_sig) * .001;
	else        
        peak_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
        peak_Integration_time_sig = (j_sig - i_sig) * .001;
    end
       
    for i_zerocrossing_rebound = i_zerocrossing:-1:1
    	if ~(splinedSTA(i_zerocrossing_rebound) < mean(trial_estim_amps)); break; end
    end
    
    plot(STA_t(1) - .001 + (i_zerocrossing_rebound * .001), splinedSTA(i_zerocrossing_rebound), 'ok', 'LineWidth', 2)
    plot(STA_t(1) - .001 + ((i_zerocrossing + 1) * .001), splinedSTA(i_zerocrossing + 1), 'ok', 'LineWidth', 2)
    
    [peak_amp_rebound, peak_location_rebound] = min(splinedSTA(i_zerocrossing_rebound:i_zerocrossing));
    
    [H_peak_rebound, ~] = ztest(peak_amp_rebound, STA_mean, STA_std, 'Alpha',alpha_value);
    
    if H_peak_rebound
        if ~ p.vstim
            peak_Integration_time = (i_zerocrossing - i_zerocrossing_rebound) * .001;
        elseif p.vstim
            trough_Integration_time = (i_zerocrossing - i_zerocrossing_rebound) * .001;
        end
        peak_location_time_rebound = abs(splinedSTA_t(1) + peak_location_rebound * .001 - .001);
        trough_location_time_rebound = - 10;
    else
        if ~ p.vstim
            peak_Integration_time = 0;
            peak_location_time_rebound = - 10;
            trough_location_time_rebound = - 10;
        elseif p.vstim
            trough_Integration_time = 0;
        end
    end
	
elseif (splinedSTA_peak_idx > splinedSTA_trough_idx)
    
    locs = find(splinedSTA(1:length(splinedSTA) / 2) < mean(trial_estim_amps));
    i_zerocrossing = locs(locs == splinedSTA_peak_idx);
    j_zerocrossing = locs(locs == splinedSTA_peak_idx);
    
    i_sig = i_zerocrossing;
    j_sig = j_zerocrossing;
    
    while (splinedSTA(i_zerocrossing) < mean(trial_estim_amps))
        i_zerocrossing = i_zerocrossing - 1;
        if (i_zerocrossing == 0)
            i_zerocrossing = 1;
            break;
        end
    end
    
    while ((splinedSTA(i_sig) < mean(trial_estim_amps)) && ztest(splinedSTA(i_sig), STA_mean, STA_std, 'Alpha',alpha_value))
        i_sig = i_sig - 1;
        if (i_sig == 0)
            i_sig = 1;
            break;
        end
    end
    
    plot(STA_t(1) - .001 + (i_zerocrossing * .001), splinedSTA(i_zerocrossing), 'ok', 'LineWidth', 2)
    hold on
    plot(STA_t(1) - .001 + (i_sig * .001), splinedSTA(i_sig), 'or', 'LineWidth', 2)
    
    while (splinedSTA(j_zerocrossing) < mean(trial_estim_amps))
        j_zerocrossing = j_zerocrossing + 1;
        if (j_zerocrossing > length(splinedSTA) / 2)
            break;
        end
        
    end
    
    while ((splinedSTA(j_sig) < mean(trial_estim_amps)) && ztest(splinedSTA(j_sig), STA_mean, STA_std, 'Alpha',alpha_value))
        j_sig = j_sig + 1;
        if (j_sig > length(splinedSTA) / 2)
            break;
        end
    end
    
    plot(STA_t(1) - .001 + (j_zerocrossing * .001), splinedSTA(j_zerocrossing), 'ok', 'LineWidth', 2)
    hold on
    plot(STA_t(1) - .001 + (j_sig * .001), splinedSTA(j_sig), 'or', 'LineWidth', 2)
    
    if ~ p.vstim
        
        peak_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
        peak_Integration_time_sig = (j_sig - i_sig) * .001;
        
    elseif p.vstim
        
        trough_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
        trough_Integration_time_sig = (j_sig - i_sig) * .001;
        
    end
    
    i_zerocrossing_rebound = i_zerocrossing;
    j_zerocrossing_rebound_end = i_zerocrossing_rebound;
    
    while (splinedSTA(i_zerocrossing_rebound) > mean(trial_estim_amps))
        i_zerocrossing_rebound = i_zerocrossing_rebound - 1;
        if (i_zerocrossing_rebound == 1)
            break;
        end
    end
    
    plot(STA_t(1) - .001 + (i_zerocrossing_rebound * .001), splinedSTA(i_zerocrossing_rebound), 'ok', 'LineWidth', 2)
    plot(STA_t(1) - .001 + ((i_zerocrossing + 1) * .001), splinedSTA(i_zerocrossing + 1), 'ok', 'LineWidth', 2)
    
    trough_amp_rebound = max(splinedSTA(i_zerocrossing_rebound:j_zerocrossing_rebound_end));
    
    [H_trough_rebound P_trough_rebound] = ztest(trough_amp_rebound, STA_mean, STA_std, 'Alpha',alpha_value);
    
    if H_trough_rebound
        
        if ~ p.vstim
            trough_Integration_time = (j_zerocrossing_rebound_end - i_zerocrossing_rebound) * .001;
        elseif p.vstim
            peak_Integration_time = (j_zerocrossing_rebound_end - i_zerocrossing_rebound) * .001;
        end
        
        trough_location_rebound = find(splinedSTA(1:length(splinedSTA) / 2) == trough_amp_rebound);
        trough_location_time_rebound = abs(splinedSTA_t(1) + trough_location_rebound * .001 - .001);
        peak_location_time_rebound = - 10;
        
    else
        
        if ~ p.vstim
            trough_Integration_time = 0;
            peak_location_time_rebound = - 10;
            trough_location_time_rebound = - 10;
        elseif p.vstim
            peak_Integration_time = 0;
        end
    end
    
    if STA_through_issig == 1
        
        locs = find(splinedSTA(1:length(splinedSTA) / 2) > mean(trial_estim_amps));
        i_zerocrossing = (locs(find(locs == splinedSTA_trough_idx)));
        j_zerocrossing = (locs(find(locs == splinedSTA_trough_idx)));
        i_sig = i_zerocrossing;
        j_sig = j_zerocrossing;
        
        while ((splinedSTA(i_sig) > mean(trial_estim_amps)) && ztest(splinedSTA(i_sig), STA_mean, STA_std, 'Alpha',alpha_value))
            i_sig = i_sig - 1;
            if (i_sig == 0)
                
                i_sig = 1;
                break;
            end
        end
        
        while ((splinedSTA(j_sig) > mean(trial_estim_amps)) && ztest(splinedSTA(j_sig), STA_mean, STA_std, 'Alpha',alpha_value))
            j_sig = j_sig + 1;
            if (j_sig > length(splinedSTA) / 2)
                break;
            end
        end
        
        trough_Integration_time_sig = (j_sig - i_sig) * .001;
        
        plot(STA_t(1) - .001 + (j_sig * .001), splinedSTA(j_sig), 'or', 'LineWidth', 2)
        hold on
        plot(STA_t(1) - .001 + (i_sig * .001), splinedSTA(i_sig), 'or', 'LineWidth', 2)
        
    else
        trough_Integration_time_sig = 0;
    end
    
elseif (splinedSTA_peak_idx == splinedSTA_trough_idx)
    peak_Integration_time = 0;
    trough_Integration_time = 0;
end

%% p.single_pulse_activation_correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~ p.vstim
    
    if p.single_pulse_activation_correction
        
        if ((((peak_location_time < .04) && (STA_peak_issig == 1)) || (((trough_location_time < .04) && (STA_through_issig == 1)))))
            disp('single pulse')
            
            correctedSTA(25) = mean(estim_meanline);
            figure;
            title('check')
            splinedSTA = spline(STA_t, correctedSTA, splinedSTA_t);
            
            plot(splinedSTA_t, splinedSTA, 'LineWidth', 2)
            ylim([plt_ylim(1) plt_ylim(2)])
            set(gcf, 'color', 'w');
            estim_meanline = (1 * mean(trial_estim_amps)) * ones(length(splinedSTA), 1);
            hold on
            plot(splinedSTA_t, estim_meanline, 'k')
            step_size = (plt_ylim(2) - plt_ylim(1)) / 10;
            zeromarker = zeros(length(plt_ylim(1):step_size:plt_ylim(2)));
            plot(zeromarker, plt_ylim(1):step_size:plt_ylim(2), 'k');
            
            %% significance of peak and trough of STA
            
            STA_mean = mean(STA(length(STA) / 2 + 1:end));
            STA_std = std(STA(length(STA) / 2 + 1:end));
            
            alpha_value = 1 - nthroot(.95, p.tKerLen);
            
            alpha_value = alpha_value / 2;
            
            STA_peak_val = min(STA(1:length(STA) / 2));
            [STA_peak_issig P_peak] = ztest(STA_peak_val, STA_mean, STA_std, 'Alpha',alpha_value);
            
            STA_trough_val = max(STA(1:length(STA) / 2));
            [STA_through_issig P_trough] = ztest(STA_trough_val, STA_mean, STA_std, 'Alpha',alpha_value);
            
            %% Location of peak and trough and their values
            
            if ~ p.vstim
                
                splinedSTA_peak_val = min(splinedSTA(1:(length(splinedSTA) / 2)));
                
            elseif p.vstim
                
                splinedSTA_peak_val = min(splinedSTA(1:(length(splinedSTA) / 2)));
                
            end
            
            % trough_val = max(splinedSTA(length(splinedSTA)/4:length(splinedSTA)/2));
            splinedSTA_trough_val = max(splinedSTA(1:length(splinedSTA) / 2));
            
            baseline_avg_right = mean(splinedSTA(1 + (length(splinedSTA) / 2):length(splinedSTA)));
            baseline_std_right = std(splinedSTA(1 + (length(splinedSTA) / 2):length(splinedSTA)));
            
            splinedSTA_trough_idx = find(splinedSTA(1:length(splinedSTA) / 2) == splinedSTA_trough_val);
            splinedSTA_peak_idx = find(splinedSTA(1:length(splinedSTA) / 2) == splinedSTA_peak_val);
            
            if ~ p.vstim
                trough_location_time = abs(splinedSTA_t(1) + splinedSTA_trough_idx * .001 - .001);
                peak_location_time = abs(splinedSTA_t(1) + splinedSTA_peak_idx * .001 - .001);
                
            elseif p.vstim
                
                trough_location_time = abs(splinedSTA_t(1) + splinedSTA_peak_idx * .001 - .001);
                peak_location_time = abs(splinedSTA_t(1) + splinedSTA_trough_idx * .001 - .001);
                
            end
            
            %% integration time calculation peak
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            if STA_peak_issig == 1
                
                locs = find(splinedSTA(1:length(splinedSTA) / 2) < mean(trial_estim_amps));
                i_zerocrossing = (locs(find(locs == splinedSTA_peak_idx)));
                j_zerocrossing = (locs(find(locs == splinedSTA_peak_idx)));
                splinedSTA_peak_idx = i_zerocrossing;
                
            else
                
                splinedSTA_peak_idx = - 1000;
                
            end
            
            if STA_through_issig == 1
                
                locs = find(splinedSTA(1:length(splinedSTA) / 2) > mean(trial_estim_amps));
                i_zerocrossing = (locs(find(locs == splinedSTA_trough_idx)));
                j_zerocrossing = (locs(find(locs == splinedSTA_trough_idx)));
                splinedSTA_trough_idx = i_zerocrossing;
                
            else
                
                splinedSTA_trough_idx = - 1000;
                
            end
            
            if (splinedSTA_peak_idx > splinedSTA_trough_idx)
                
                locs = find(splinedSTA(1:length(splinedSTA) / 2) < mean(trial_estim_amps));
                i_zerocrossing = (locs(find(locs == splinedSTA_peak_idx)));
                j_zerocrossing = (locs(find(locs == splinedSTA_peak_idx)));
                
                i_sig = i_zerocrossing;
                j_sig = j_zerocrossing;
                
                while (splinedSTA(i_zerocrossing) < mean(trial_estim_amps))
                    
                    i_zerocrossing = i_zerocrossing - 1;
                    if (i_zerocrossing == 0)
                        i_zerocrossing = 1;
                        break;
                    end
                    
                end
                
                while ((splinedSTA(i_sig) < mean(trial_estim_amps)) && ztest(splinedSTA(i_sig), STA_mean, STA_std, 'Alpha',alpha_value))
                    
                    i_sig = i_sig - 1;
                    if (i_sig == 0)
                        i_sig = 1;
                        break;
                    end
                    
                end
                
                plot(STA_t(1) - .001 + (i_zerocrossing * .001), splinedSTA(i_zerocrossing), 'ok', 'LineWidth', 2)
                hold on
                plot(STA_t(1) - .001 + (i_sig * .001), splinedSTA(i_sig), 'or', 'LineWidth', 2)
                
                while (splinedSTA(j_zerocrossing) < mean(trial_estim_amps))
                    j_zerocrossing = j_zerocrossing + 1;
                    if (j_zerocrossing > length(splinedSTA) / 2)
                        break;
                    end
                    
                end
                
                while ((splinedSTA(j_sig) < mean(trial_estim_amps)) && ztest(splinedSTA(j_sig), STA_mean, STA_std, 'Alpha',alpha_value))
                    
                    j_sig = j_sig + 1;
                    if (j_sig > length(splinedSTA) / 2)
                        break;
                    end
                    
                end
                
                plot(STA_t(1) - .001 + (j_zerocrossing * .001), splinedSTA(j_zerocrossing), 'ok', 'LineWidth', 2)
                hold on
                plot(STA_t(1) - .001 + (j_sig * .001), splinedSTA(j_sig), 'or', 'LineWidth', 2)
                
                if ~ p.vstim
                    
                    peak_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
                    peak_Integration_time_sig = (j_sig - i_sig) * .001;
                    
                elseif p.vstim
                    
                    trough_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
                    trough_Integration_time_sig = (j_sig - i_sig) * .001;
                    
                end
                
                i_zerocrossing_rebound = i_zerocrossing;
                j_zerocrossing_rebound_end = i_zerocrossing_rebound;
                
                while (splinedSTA(i_zerocrossing_rebound) > mean(trial_estim_amps))
                    i_zerocrossing_rebound = i_zerocrossing_rebound - 1;
                    if (i_zerocrossing_rebound == 1)
                        break;
                    end
                end
                
                plot(STA_t(1) - .001 + (i_zerocrossing_rebound * .001), splinedSTA(i_zerocrossing_rebound), 'ok', 'LineWidth', 2)
                plot(STA_t(1) - .001 + ((i_zerocrossing + 1) * .001), splinedSTA(i_zerocrossing + 1), 'ok', 'LineWidth', 2)
                
                trough_amp_rebound = max(splinedSTA(i_zerocrossing_rebound:j_zerocrossing_rebound_end));
                
                [H_trough_rebound P_trough_rebound] = ztest(trough_amp_rebound, STA_mean, STA_std, 'Alpha',alpha_value);
                
                if H_trough_rebound
                    
                    if ~ p.vstim
                        trough_Integration_time = (j_zerocrossing_rebound_end - i_zerocrossing_rebound) * .001;
                        
                    elseif p.vstim
                        
                        peak_Integration_time = (j_zerocrossing_rebound_end - i_zerocrossing_rebound) * .001;
                        
                    end
                    
                    trough_location_rebound = find(splinedSTA(1:length(splinedSTA) / 2) == trough_amp_rebound);
                    trough_location_time_rebound = abs(splinedSTA_t(1) + trough_location_rebound * .001 - .001);
                    peak_location_time_rebound = - 10;
                    
                else
                    
                    if ~ p.vstim
                        
                        trough_Integration_time = 0;
                        peak_location_time_rebound = - 10;
                        trough_location_time_rebound = - 10;
                        
                    elseif p.vstim
                        
                        peak_Integration_time = 0;
                        
                    end
                    
                end
                
                if STA_through_issig == 1
                    
                    locs = find(splinedSTA(1:length(splinedSTA) / 2) > mean(trial_estim_amps));
                    i_zerocrossing = (locs(find(locs == splinedSTA_trough_idx)));
                    j_zerocrossing = (locs(find(locs == splinedSTA_trough_idx)));
                    i_sig = i_zerocrossing;
                    j_sig = j_zerocrossing;
                    
                    while ((splinedSTA(i_sig) > mean(trial_estim_amps)) && ztest(splinedSTA(i_sig), STA_mean, STA_std, 'Alpha',alpha_value))
                        
                        i_sig = i_sig - 1;
                        if (i_sig == 0)
                            
                            i_sig = 1;
                            break;
                        end
                        
                    end
                    
                    while ((splinedSTA(j_sig) > mean(trial_estim_amps)) && ztest(splinedSTA(j_sig), STA_mean, STA_std, 'Alpha',alpha_value))
                        
                        j_sig = j_sig + 1;
                        if (j_sig > length(splinedSTA) / 2)
                            break;
                        end
                        
                    end
                    
                    trough_Integration_time_sig = (j_sig - i_sig) * .001;
                    plot(STA_t(1) - .001 + (j_sig * .001), splinedSTA(j_sig), 'or', 'LineWidth', 2)
                    hold on
                    plot(STA_t(1) - .001 + (i_sig * .001), splinedSTA(i_sig), 'or', 'LineWidth', 2)
                    
                else
                    
                    trough_Integration_time_sig = 0;
                    
                end
                
            elseif (splinedSTA_peak_idx < splinedSTA_trough_idx)
                
                locs = find(splinedSTA(1:length(splinedSTA) / 2) > mean(trial_estim_amps));
                i_zerocrossing = (locs(find(locs == splinedSTA_trough_idx)));
                j_zerocrossing = (locs(find(locs == splinedSTA_trough_idx)));
                i_sig = i_zerocrossing;
                j_sig = j_zerocrossing;
                
                while (splinedSTA(i_zerocrossing) > mean(trial_estim_amps))
                    
                    i_zerocrossing = i_zerocrossing - 1;
                    
                    if (i_zerocrossing == 0)
                        i_zerocrossing = 1;
                        break;
                        
                    end
                    
                end
                
                while ((splinedSTA(i_sig) > mean(trial_estim_amps)) && ztest(splinedSTA(i_sig), STA_mean, STA_std, 'Alpha',alpha_value))
                    
                    i_sig = i_sig - 1;
                    if (i_sig == 0)
                        i_sig = 1;
                        break;
                    end
                    
                end
                
                plot(STA_t(1) - .001 + (i_zerocrossing * .001), splinedSTA(i_zerocrossing), 'ok', 'LineWidth', 2)
                hold on
                plot(STA_t(1) - .001 + (i_sig * .001), splinedSTA(i_sig), 'or', 'LineWidth', 2)
                
                while (splinedSTA(j_zerocrossing) > mean(trial_estim_amps))
                    j_zerocrossing = j_zerocrossing + 1;
                    if (j_zerocrossing > length(splinedSTA) / 2)
                        break;
                    end
                end
                
                while ((splinedSTA(j_sig) > mean(trial_estim_amps)) && ztest(splinedSTA(j_sig), STA_mean, STA_std, 'Alpha',alpha_value))
                    
                    j_sig = j_sig + 1;
                    if (j_sig > length(splinedSTA) / 2)
                        break;
                    end
                end
                
                plot(STA_t(1) - .001 + (j_zerocrossing * .001), splinedSTA(j_zerocrossing), 'ok', 'LineWidth', 2)
                hold on
                plot(STA_t(1) - .001 + (j_sig * .001), splinedSTA(j_sig), 'or', 'LineWidth', 2)
                
                if ~ p.vstim
                    
                    trough_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
                    trough_Integration_time_sig = (j_sig - i_sig) * .001;
                    
                elseif p.vstim
                    
                    peak_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
                    peak_Integration_time_sig = (j_sig - i_sig) * .001;
                    
                end
                
                i_zerocrossing_rebound = i_zerocrossing;
                j_zerocrossing_rebound_end = i_zerocrossing_rebound;
                
                while (splinedSTA(i_zerocrossing_rebound) < mean(trial_estim_amps))
                    i_zerocrossing_rebound = i_zerocrossing_rebound - 1;
                    if (i_zerocrossing_rebound == 1)
                        break;
                    end
                end
                
                plot(STA_t(1) - .001 + (i_zerocrossing_rebound * .001), splinedSTA(i_zerocrossing_rebound), 'ok', 'LineWidth', 2)
                plot(STA_t(1) - .001 + ((i_zerocrossing + 1) * .001), splinedSTA(i_zerocrossing + 1), 'ok', 'LineWidth', 2)
                
                peak_amp_rebound = min(splinedSTA(i_zerocrossing_rebound:j_zerocrossing_rebound_end));
                
                [H_peak_rebound P_peak_rebound] = ztest(peak_amp_rebound, STA_mean, STA_std, 'Alpha',alpha_value);
                
                if H_peak_rebound
                    
                    if ~ p.vstim
                        peak_Integration_time = (j_zerocrossing_rebound_end - i_zerocrossing_rebound) * .001;
                        
                    elseif p.vstim
                        
                        trough_Integration_time = (j_zerocrossing_rebound_end - i_zerocrossing_rebound) * .001;
                        
                    end
                    
                    peak_location_rebound = find(splinedSTA(1:length(splinedSTA) / 2) == peak_amp_rebound);
                    peak_location_time_rebound = abs(splinedSTA_t(1) + peak_location_rebound * .001 - .001);
                    trough_location_time_rebound = - 10;
                    
                else
                    
                    if ~ p.vstim
                        
                        peak_Integration_time = 0;
                        peak_location_time_rebound = - 10;
                        trough_location_time_rebound = - 10;
                        
                    elseif p.vstim
                        
                        trough_Integration_time = 0;
                        
                    end
                    
                end
                
                if STA_peak_issig == 1
                    
                    locs = find(splinedSTA(1:length(splinedSTA) / 2) < mean(trial_estim_amps));
                    i_zerocrossing = (locs(find(locs == splinedSTA_peak_idx)));
                    j_zerocrossing = (locs(find(locs == splinedSTA_peak_idx)));
                    i_sig = i_zerocrossing;
                    j_sig = j_zerocrossing;
                    
                    while ((splinedSTA(i_sig) < mean(trial_estim_amps)) && ztest(splinedSTA(i_sig), STA_mean, STA_std, 'Alpha',alpha_value))
                        
                        i_sig = i_sig - 1;
                        if (i_sig == 0)
                            
                            i_sig = 1;
                            break;
                        end
                        
                    end
                    
                    while ((splinedSTA(j_sig) < mean(trial_estim_amps)) && ztest(splinedSTA(j_sig), STA_mean, STA_std, 'Alpha',alpha_value))
                        
                        j_sig = j_sig + 1;
                        if (j_sig > length(splinedSTA) / 2)
                            break;
                        end
                        
                    end
                    
                    peak_Integration_time_sig = (j_sig - i_sig) * .001;
                    plot(STA_t(1) - .001 + (j_sig * .001), splinedSTA(j_sig), 'or', 'LineWidth', 2)
                    hold on
                    plot(STA_t(1) - .001 + (i_sig * .001), splinedSTA(i_sig), 'or', 'LineWidth', 2)
                    
                else
                    
                    peak_Integration_time_sig = 0;
                    
                end
                
            elseif (splinedSTA_peak_idx == splinedSTA_trough_idx)
                
                peak_Integration_time = 0;
                trough_Integration_time = 0;
                
            end
            
        end
        
    end
    
end

%%

if ~ p.vstim
    
    text(- .8, - 350, 'Peak Location is ')
    text(- .8, - 385, num2str(peak_location_time_rebound))
    text(- .8, - 415, num2str(peak_location_time))
    
    text(- .8, - 450, 'Peak Integration Time is ')
    text(- .8, - 485, num2str(peak_Integration_time))
    text(- .8, - 515, num2str(peak_Integration_time_sig))
    
    text(- .8, - 550, 'Trough Location is ')
    text(- .8, - 585, num2str(trough_location_time_rebound))
    text(- .8, - 615, num2str(trough_location_time))
    
    text(- .8, - 650, 'Trough Integration Time is ')
    text(- .8, - 685, num2str(trough_Integration_time))
    text(- .8, - 715, num2str(trough_Integration_time_sig))
    
elseif p.vstim
    
    text(.4, .9, 'Peak Location is ')
    text(.4, .8, num2str(peak_location_time))
    
    text(.4, .7, 'Peak Integration Time is ')
    text(.4, .6, num2str(peak_Integration_time))
    
    text(.4, .5, 'Trough Location is ')
    text(.4, .4, num2str(trough_location_time))
    
    text(.4, .3, 'Trough Integration Time is ')
    text(.4, .2, num2str(trough_Integration_time))
    
end

if (trough_Integration_time > 0 && peak_Integration_time > 0)
    if peak_location_time < trough_location_time
        
        cell_type = 'ON';
        
    elseif peak_location_time > trough_location_time
        
        cell_type = 'OFF';
        
    end
    
elseif (trough_Integration_time == 0)
    
    cell_type = 'ON';
    
elseif (peak_Integration_time == 0)
    cell_type = 'OFF';
    
end

title(cell_type)

peak_location_time
peak_Integration_time_zero_crossing = peak_Integration_time;
trough_location_time
trough_Integration_time_zero_crossing = trough_Integration_time;

saveas(gcf, [p.work_dir, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_trial), '-', 'cubic splined', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), '.fig'], 'fig');
set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
saveas(gcf, [p.work_dir, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_trial), '-', 'cubic splined', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), '.jpeg'], 'jpeg');

%% p.weighted_burst, p.singleton_spikes
date1 = p.year;
date1(end) = [];

if ~ p.cardinal_STA_Only_Burst
    
    p.work_dir = fullfile(p.work_dir, 'Population_Analysis_Alternate\');
    if ~ p.vstim
        file_loc = strcat(p.work_dir, p.cell_id, date1, 'single_pulse_analysis = ', num2str(p.single_pulse_activation_correction));
    elseif p.vstim
        file_loc = strcat(p.work_dir, p.cell_id, date1);
    end
    
else
    
    p.work_dir = fullfile(p.work_dir, 'Population_Analysis_Alternate_BSTA\');
    
    file_loc = strcat(p.work_dir, p.cell_id, date1, ' p.weighted_burst = ', num2str(p.weighted_burst), ' p.singleton_spikes = ', num2str(p.singleton_spikes), 'single_pulse_analysis = ', num2str(p.single_pulse_activation_correction));
    
end

if p.vstim
    
    temp_store = splinedSTA_peak_val;
    splinedSTA_peak_val = splinedSTA_trough_val;
    splinedSTA_trough_val = temp_store;
    
end

if ~ exist(p.work_dir, 'dir'), mkdir(p.work_dir); end
save(file_loc, 'peak_location_time', 'trough_location_time', 'cell_type', 'STA_peak_val', 'STA_trough_val', 'peak_Integration_time_zero_crossing', 'trough_Integration_time_zero_crossing', 'peak_Integration_time_sig', 'trough_Integration_time_sig', 'peak_location_time_rebound', 'trough_location_time_rebound')

end