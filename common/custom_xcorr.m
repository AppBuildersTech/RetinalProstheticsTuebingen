function [out_xcorr] = custom_xcorr(Xin,kernel,xcorr_type)
% computed the cross corelation between a 1D signal and a kernel
% Xin: 1 x T signal
% kernel: 1 x Kw kernel
if nargin<3
    xcorr_type = 'valid';
end

Kw = size(kernel,1);

if Kw<2;error('The kernel shape should be 1 x Kw!');end

if strcmp(xcorr_type,'full')
    Xin = vertcat(zeros(Kw-1,1), Xin,zeros(Kw-1,1));
end

T = size(Xin,1);

out_xcorr = zeros(T-Kw+1,1);

for xIdx = 1:T-Kw+1
    out_xcorr(xIdx,1) = Xin(xIdx:(xIdx+Kw-1),1)' * kernel;
end

end

