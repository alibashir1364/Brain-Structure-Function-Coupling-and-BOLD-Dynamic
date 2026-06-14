function RLFP = compute_RLFP(X, TR)
% Relative Low-Frequency Power
[N, T] = size(X);
Fs = 1/TR;

freqs = (0:T-1)*(Fs/T);
low_idx = freqs < 0.14;

RLFP = zeros(N,1);

for i = 1:N
    x = detrend(X(i,:));
    Pxx = abs(fft(x)).^2;
    RLFP(i) = sum(Pxx(low_idx)) / sum(Pxx);
end
end
