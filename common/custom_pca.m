function [Um, variance_perserved]=custom_pca(data, nPCs)

Cov_data = cov(data);

[U, S, ~] = svd(Cov_data);

eigvals = diag(S); % these are square of variances in the aligned to the maximum variations space 
N = size(data,2);
sumeigvals = sum(eigvals);
variance_sum_upto_m = zeros(1,N);
for m=1:N
    variance_sum_upto_m(m) = sum(abs(eigvals(1:m)));
end
variance_sum_upto_m = variance_sum_upto_m/sumeigvals;
variance_perserved = variance_sum_upto_m(nPCs);

Um = U(:,1:nPCs);
      
end
