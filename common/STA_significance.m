function out_p = STA_significance(STA_LS, STA_RS, estim_mean, exp_ps)
    % Signifance Algorithem
    % 	1 - Find the min/max deflections in the magnitude of the splined STA
    % 	2 - Check for significance of the deflections
    %   3 - Find the crossings of the STA with the baseline after/before the previous computed deflectios
    % 	4 - Assign the one that comes closer to time zero (closer to half the STA length) to D1, and the other to D2

    % 	1 - Find the min/max deflections in the magnitude of the splined STA
    [~, maxD_idx] = max(sqrt(STA_LS.^2));
    [~, minD_idx] = min(sqrt(STA_LS.^2));

    maxD_val = STA_LS(maxD_idx);
    minD_val = STA_LS(minD_idx);

    % 	2 - Check for significance of the deflections
    alpha_val = (1 - nthroot(.95, exp_ps.tKerLen))/2;

    STA_RS_mean = mean(STA_RS);%right side mean
    STA_RS_std = std(STA_RS);%right side std

    [maxD_issig, ~] = ztest(maxD_val, STA_RS_mean, STA_RS_std, 'Alpha', alpha_val); % h = 1 means an outlier, and we want an outlier
    [minD_issig, ~] = ztest(minD_val, STA_RS_mean, STA_RS_std, 'Alpha', alpha_val);

    %   3 - Find the crossings of the STA with the baseline after/before the previous computed deflectios
    [maxD_zcross_p1,maxD_zcross_p2] = find_crossing(STA_LS, maxD_idx, estim_mean);
    [minD_zcross_p1,minD_zcross_p2] = find_crossing(STA_LS, minD_idx, estim_mean);

    [maxD_finsig_p1,maxD_finsig_p2] = find_first_insignificance(STA_LS, maxD_idx, STA_RS_mean, STA_RS_std, alpha_val);
    [minD_finsig_p1,minD_finsig_p2] = find_first_insignificance(STA_LS, minD_idx, STA_RS_mean, STA_RS_std, alpha_val);

    % 	4 - Assign the one that comes closer to time zero (closer to half the STA length) to D1, and the other to D2
    if maxD_idx > minD_idx
        out_p.D1_idx = maxD_idx;
        out_p.D1_val = maxD_val;
        out_p.D1_issig = maxD_issig;
        out_p.D1_cross_ids = [maxD_zcross_p1,maxD_zcross_p2];
        out_p.D1_finsig_ids = [maxD_finsig_p1,maxD_finsig_p2];

        out_p.D2_idx = minD_idx;
        out_p.D2_val = minD_val;
        out_p.D2_issig = minD_issig;
        out_p.D2_cross_ids = [minD_zcross_p1,minD_zcross_p2];
        out_p.D2_finsig_ids = [minD_finsig_p1,minD_finsig_p2];
    else
        out_p.D1_idx = minD_idx;
        out_p.D1_val = minD_val;
        out_p.D1_issig = minD_issig;
        out_p.D1_cross_ids = [minD_zcross_p1,minD_zcross_p2];
        out_p.D1_finsig_ids = [minD_finsig_p1,minD_finsig_p2];

        out_p.D2_idx = maxD_idx;
        out_p.D2_val = maxD_val;
        out_p.D2_issig = maxD_issig;
        out_p.D2_cross_ids = [maxD_zcross_p1,maxD_zcross_p2];
        out_p.D2_finsig_ids = [maxD_finsig_p1,maxD_finsig_p2];
    end
    %trough_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
end

function [p1,p2] = find_crossing(signal_vals, peak_idx, baseline)
    p1 = NaN; p2 = NaN;
    for sIdx = peak_idx:-1:1
        if signal_vals(peak_idx) > baseline
            if signal_vals(sIdx)<=baseline
                p1 = sIdx;
                break;
            end
        end
        if signal_vals(peak_idx) <= baseline
            if signal_vals(sIdx)>= baseline
                p1 = sIdx;
                break;
            end
        end
    end
    for sIdx = peak_idx:length(signal_vals)
        if signal_vals(peak_idx) > baseline
            if signal_vals(sIdx)<=baseline
                p2 = sIdx;
                break;
            end
        end
        if signal_vals(peak_idx) <= baseline
            if signal_vals(sIdx)>= baseline
                p2 = sIdx;
                break;
            end
        end
    end
end
function [p1,p2] = find_first_insignificance(signal_vals, peak_idx, z_mean, z_std, alpha_val)
    % finds the first insignifanct points of the deflection after the peak value
    p1 = NaN; p2 = NaN;
    for sIdx = peak_idx:-1:1
        if ~ztest(signal_vals(sIdx), z_mean, z_std, 'Alpha', alpha_val)
            p1 = sIdx;
            break;
        end
    end
    for sIdx = peak_idx:length(signal_vals)
        if ~ztest(signal_vals(sIdx), z_mean, z_std, 'Alpha', alpha_val)
            p2 = sIdx;
            break;
        end
    end
end