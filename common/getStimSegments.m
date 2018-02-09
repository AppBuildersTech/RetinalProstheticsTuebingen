function [raw_stim_ensem, sp_assoc_ensem] = getStimSegments(estim,estim_ts,estim_spts, Kw)
    raw_stim_ensem = zeros(T-Kw+1,Kw);
    sp_assoc_ensem = [];

    for xIdx = 1:T-Kw+1
        raw_stim_ensem(xIdx,:) = estim(xIdx:(xIdx+Kw-1),1);
        window_end_t = estim_ts(xIdx+Kw-1,1);
        sp_ids = estim_spts>=window_end_t & estim_spts<=window_end_t+speriod;
        if sum(sp_ids)>0
            sp_assoc_ensem = vertcat(sp_assoc_ensem, raw_stim_ensem(xIdx,:));
            % number of sp_assoc_estim_excerpts might be smaller than the spike counts. that would happen if the few first spikes actually happened within the first window
        end
    end
end