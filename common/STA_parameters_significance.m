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

%% cubic splined sta


%% significance of peak and trough of STA

baseline_sta = mean(STA(length(STA) / 2 + 1:end));
std_baselines = std(STA(length(STA) / 2 + 1:end));

alphaville = 1 - nthroot(.95, p.tKerLen);

alphaville = alphaville / 2;

peak_sta = min(STA(1:length(STA) / 2));
[H_peak P_peak] = ztest(peak_sta, baseline_sta, std_baselines, alphaville);

trough_sta = max(STA(1:length(STA) / 2));
[H_trough P_trough] = ztest(trough_sta, baseline_sta, std_baselines, alphaville);

%% Location of peak and trough and their values

if ~ p.vstim
    
    peak_T = min(splinedSTA(1:(length(splinedSTA) / 2)));
    
elseif p.vstim
    
    peak_T = min(splinedSTA(1:(length(splinedSTA) / 2)));
    
end

% trough_T = max(splinedSTA(length(splinedSTA)/4:length(splinedSTA)/2));
trough_T = max(splinedSTA(1:length(splinedSTA) / 2));

baseline_avg_right = mean(splinedSTA(1 + (length(splinedSTA) / 2):length(splinedSTA)));
baseline_std_right = std(splinedSTA(1 + (length(splinedSTA) / 2):length(splinedSTA)));

trough_location = find(splinedSTA(1:length(splinedSTA) / 2) == trough_T);
peak_location = find(splinedSTA(1:length(splinedSTA) / 2) == peak_T);

if ~ p.vstim
    trough_location_time = abs(splinedSTA_t(1) + trough_location * .001 - .001);
    peak_location_time = abs(splinedSTA_t(1) + peak_location * .001 - .001);
    
elseif p.vstim
    
    trough_location_time = abs(splinedSTA_t(1) + peak_location * .001 - .001);
    peak_location_time = abs(splinedSTA_t(1) + trough_location * .001 - .001);
    
end

%% integration time calculation peak

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if H_peak == 1
    
    locs = find(splinedSTA(1:length(splinedSTA) / 2) < mean(trial_estim_amps));
    i_zerocrossing = (locs(find(locs == peak_location)));
    j_zerocrossing = (locs(find(locs == peak_location)));
    i_peak = i_zerocrossing;
    
else
    
    i_peak = - 1000;
    
end

if H_trough == 1
    
    locs = find(splinedSTA(1:length(splinedSTA) / 2) > mean(trial_estim_amps));
    i_zerocrossing = (locs(find(locs == trough_location)));
    j_zerocrossing = (locs(find(locs == trough_location)));
    i_trough = i_zerocrossing;
    
else
    
    i_trough = - 1000;
    
end

if (i_peak > i_trough)
    
    locs = find(splinedSTA(1:length(splinedSTA) / 2) < mean(trial_estim_amps));
    i_zerocrossing = (locs(find(locs == peak_location)));
    j_zerocrossing = (locs(find(locs == peak_location)));
    
    i_sig = i_zerocrossing;
    j_sig = j_zerocrossing;
    
    while (splinedSTA(i_zerocrossing) < mean(trial_estim_amps))
        
        i_zerocrossing = i_zerocrossing - 1;
        if (i_zerocrossing == 0)
            i_zerocrossing = 1;
            break;
        end
        
    end
    
    while ((splinedSTA(i_sig) < mean(trial_estim_amps)) && ztest(splinedSTA(i_sig), baseline_sta, std_baselines, alphaville))
        
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
    
    while ((splinedSTA(j_sig) < mean(trial_estim_amps)) && ztest(splinedSTA(j_sig), baseline_sta, std_baselines, alphaville))
        
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
    
    [H_trough_rebound P_trough_rebound] = ztest(trough_amp_rebound, baseline_sta, std_baselines, alphaville);
    
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
    
    if H_trough == 1
        
        locs = find(splinedSTA(1:length(splinedSTA) / 2) > mean(trial_estim_amps));
        i_zerocrossing = (locs(find(locs == trough_location)));
        j_zerocrossing = (locs(find(locs == trough_location)));
        i_sig = i_zerocrossing;
        j_sig = j_zerocrossing;
        
        while ((splinedSTA(i_sig) > mean(trial_estim_amps)) && ztest(splinedSTA(i_sig), baseline_sta, std_baselines, alphaville))
            
            i_sig = i_sig - 1;
            if (i_sig == 0)
                
                i_sig = 1;
                break;
            end
            
        end
        
        while ((splinedSTA(j_sig) > mean(trial_estim_amps)) && ztest(splinedSTA(j_sig), baseline_sta, std_baselines, alphaville))
            
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
    
elseif (i_peak < i_trough)
    
    locs = find(splinedSTA(1:length(splinedSTA) / 2) > mean(trial_estim_amps));
    i_zerocrossing = (locs(find(locs == trough_location)));
    j_zerocrossing = (locs(find(locs == trough_location)));
    i_sig = i_zerocrossing;
    j_sig = j_zerocrossing;
    
    while (splinedSTA(i_zerocrossing) > mean(trial_estim_amps))
        
        i_zerocrossing = i_zerocrossing - 1;
        
        if (i_zerocrossing == 0)
            i_zerocrossing = 1;
            break;
            
        end
        
    end
    
    while ((splinedSTA(i_sig) > mean(trial_estim_amps)) && ztest(splinedSTA(i_sig), baseline_sta, std_baselines, alphaville))
        
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
    
    while ((splinedSTA(j_sig) > mean(trial_estim_amps)) && ztest(splinedSTA(j_sig), baseline_sta, std_baselines, alphaville))
        
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
    
    [H_peak_rebound P_peak_rebound] = ztest(peak_amp_rebound, baseline_sta, std_baselines, alphaville);
    
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
    
    if H_peak == 1
        
        locs = find(splinedSTA(1:length(splinedSTA) / 2) < mean(trial_estim_amps));
        i_zerocrossing = (locs(find(locs == peak_location)));
        j_zerocrossing = (locs(find(locs == peak_location)));
        i_sig = i_zerocrossing;
        j_sig = j_zerocrossing;
        
        while ((splinedSTA(i_sig) < mean(trial_estim_amps)) && ztest(splinedSTA(i_sig), baseline_sta, std_baselines, alphaville))
            
            i_sig = i_sig - 1;
            if (i_sig == 0)
                
                i_sig = 1;
                break;
            end
            
        end
        
        while ((splinedSTA(j_sig) < mean(trial_estim_amps)) && ztest(splinedSTA(j_sig), baseline_sta, std_baselines, alphaville))
            
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
    
elseif (i_peak == i_trough)
    
    peak_Integration_time = 0;
    trough_Integration_time = 0;
    
end

%% p.single_pulse_activation_correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~ p.vstim
    
    if p.single_pulse_activation_correction
        
        if ((((peak_location_time < .04) && (H_peak == 1)) || (((trough_location_time < .04) && (H_trough == 1)))))
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
            
            baseline_sta = mean(STA(length(STA) / 2 + 1:end));
            std_baselines = std(STA(length(STA) / 2 + 1:end));
            
            alphaville = 1 - nthroot(.95, p.tKerLen);
            
            alphaville = alphaville / 2;
            
            peak_sta = min(STA(1:length(STA) / 2));
            [H_peak P_peak] = ztest(peak_sta, baseline_sta, std_baselines, alphaville);
            
            trough_sta = max(STA(1:length(STA) / 2));
            [H_trough P_trough] = ztest(trough_sta, baseline_sta, std_baselines, alphaville);
            
            %% Location of peak and trough and their values
            
            if ~ p.vstim
                
                peak_T = min(splinedSTA(1:(length(splinedSTA) / 2)));
                
            elseif p.vstim
                
                peak_T = min(splinedSTA(1:(length(splinedSTA) / 2)));
                
            end
            
            % trough_T = max(splinedSTA(length(splinedSTA)/4:length(splinedSTA)/2));
            trough_T = max(splinedSTA(1:length(splinedSTA) / 2));
            
            baseline_avg_right = mean(splinedSTA(1 + (length(splinedSTA) / 2):length(splinedSTA)));
            baseline_std_right = std(splinedSTA(1 + (length(splinedSTA) / 2):length(splinedSTA)));
            
            trough_location = find(splinedSTA(1:length(splinedSTA) / 2) == trough_T);
            peak_location = find(splinedSTA(1:length(splinedSTA) / 2) == peak_T);
            
            if ~ p.vstim
                trough_location_time = abs(splinedSTA_t(1) + trough_location * .001 - .001);
                peak_location_time = abs(splinedSTA_t(1) + peak_location * .001 - .001);
                
            elseif p.vstim
                
                trough_location_time = abs(splinedSTA_t(1) + peak_location * .001 - .001);
                peak_location_time = abs(splinedSTA_t(1) + trough_location * .001 - .001);
                
            end
            
            %% integration time calculation peak
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            if H_peak == 1
                
                locs = find(splinedSTA(1:length(splinedSTA) / 2) < mean(trial_estim_amps));
                i_zerocrossing = (locs(find(locs == peak_location)));
                j_zerocrossing = (locs(find(locs == peak_location)));
                i_peak = i_zerocrossing;
                
            else
                
                i_peak = - 1000;
                
            end
            
            if H_trough == 1
                
                locs = find(splinedSTA(1:length(splinedSTA) / 2) > mean(trial_estim_amps));
                i_zerocrossing = (locs(find(locs == trough_location)));
                j_zerocrossing = (locs(find(locs == trough_location)));
                i_trough = i_zerocrossing;
                
            else
                
                i_trough = - 1000;
                
            end
            
            if (i_peak > i_trough)
                
                locs = find(splinedSTA(1:length(splinedSTA) / 2) < mean(trial_estim_amps));
                i_zerocrossing = (locs(find(locs == peak_location)));
                j_zerocrossing = (locs(find(locs == peak_location)));
                
                i_sig = i_zerocrossing;
                j_sig = j_zerocrossing;
                
                while (splinedSTA(i_zerocrossing) < mean(trial_estim_amps))
                    
                    i_zerocrossing = i_zerocrossing - 1;
                    if (i_zerocrossing == 0)
                        i_zerocrossing = 1;
                        break;
                    end
                    
                end
                
                while ((splinedSTA(i_sig) < mean(trial_estim_amps)) && ztest(splinedSTA(i_sig), baseline_sta, std_baselines, alphaville))
                    
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
                
                while ((splinedSTA(j_sig) < mean(trial_estim_amps)) && ztest(splinedSTA(j_sig), baseline_sta, std_baselines, alphaville))
                    
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
                
                [H_trough_rebound P_trough_rebound] = ztest(trough_amp_rebound, baseline_sta, std_baselines, alphaville);
                
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
                
                if H_trough == 1
                    
                    locs = find(splinedSTA(1:length(splinedSTA) / 2) > mean(trial_estim_amps));
                    i_zerocrossing = (locs(find(locs == trough_location)));
                    j_zerocrossing = (locs(find(locs == trough_location)));
                    i_sig = i_zerocrossing;
                    j_sig = j_zerocrossing;
                    
                    while ((splinedSTA(i_sig) > mean(trial_estim_amps)) && ztest(splinedSTA(i_sig), baseline_sta, std_baselines, alphaville))
                        
                        i_sig = i_sig - 1;
                        if (i_sig == 0)
                            
                            i_sig = 1;
                            break;
                        end
                        
                    end
                    
                    while ((splinedSTA(j_sig) > mean(trial_estim_amps)) && ztest(splinedSTA(j_sig), baseline_sta, std_baselines, alphaville))
                        
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
                
            elseif (i_peak < i_trough)
                
                locs = find(splinedSTA(1:length(splinedSTA) / 2) > mean(trial_estim_amps));
                i_zerocrossing = (locs(find(locs == trough_location)));
                j_zerocrossing = (locs(find(locs == trough_location)));
                i_sig = i_zerocrossing;
                j_sig = j_zerocrossing;
                
                while (splinedSTA(i_zerocrossing) > mean(trial_estim_amps))
                    
                    i_zerocrossing = i_zerocrossing - 1;
                    
                    if (i_zerocrossing == 0)
                        i_zerocrossing = 1;
                        break;
                        
                    end
                    
                end
                
                while ((splinedSTA(i_sig) > mean(trial_estim_amps)) && ztest(splinedSTA(i_sig), baseline_sta, std_baselines, alphaville))
                    
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
                
                while ((splinedSTA(j_sig) > mean(trial_estim_amps)) && ztest(splinedSTA(j_sig), baseline_sta, std_baselines, alphaville))
                    
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
                
                [H_peak_rebound P_peak_rebound] = ztest(peak_amp_rebound, baseline_sta, std_baselines, alphaville);
                
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
                
                if H_peak == 1
                    
                    locs = find(splinedSTA(1:length(splinedSTA) / 2) < mean(trial_estim_amps));
                    i_zerocrossing = (locs(find(locs == peak_location)));
                    j_zerocrossing = (locs(find(locs == peak_location)));
                    i_sig = i_zerocrossing;
                    j_sig = j_zerocrossing;
                    
                    while ((splinedSTA(i_sig) < mean(trial_estim_amps)) && ztest(splinedSTA(i_sig), baseline_sta, std_baselines, alphaville))
                        
                        i_sig = i_sig - 1;
                        if (i_sig == 0)
                            
                            i_sig = 1;
                            break;
                        end
                        
                    end
                    
                    while ((splinedSTA(j_sig) < mean(trial_estim_amps)) && ztest(splinedSTA(j_sig), baseline_sta, std_baselines, alphaville))
                        
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
                
            elseif (i_peak == i_trough)
                
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
    
    temp_store = peak_T;
    peak_T = trough_T;
    trough_T = temp_store;
    
end

if ~ exist(p.work_dir, 'dir'), mkdir(p.work_dir); end
save(file_loc, 'peak_location_time', 'trough_location_time', 'cell_type', 'peak_T', 'trough_T', 'baseline_avg_right', 'baseline_std_right', 'peak_Integration_time_zero_crossing', 'trough_Integration_time_zero_crossing', 'peak_Integration_time_sig', 'trough_Integration_time_sig', 'peak_location_time_rebound', 'trough_location_time_rebound')

end