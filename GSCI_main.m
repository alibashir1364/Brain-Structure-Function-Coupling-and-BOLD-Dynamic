% One-sample test vs zero
[p_low,~,stats_low]   = signrank(GSCI.GSCI_low);
[p_high,~,stats_high] = signrank(GSCI.GSCI_high);

% Paired comparison
[p_pair,~,stats_pair] = signrank(GSCI.GSCI_low, GSCI.GSCI_high);

fprintf('GSCI_low vs 0:  p = %.3e\n', p_low)
fprintf('GSCI_high vs 0: p = %.3e\n', p_high)
fprintf('Low vs High:    p = %.3e\n', p_pair)
