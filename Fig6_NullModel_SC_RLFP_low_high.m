function Fig6_NullModel_SC_RLFP_low_high(SC, RLFP_low_mean, RLFP_high_mean, roi_volume, Nperm)
% FIGURE 6 — Null Model Validation (Hypothesis H5)
% Tests whether SC–timescale coupling depends on empirical connectome
% for both coupled (low-pass) and decoupled (high-pass) signals
%
% Inputs:
%   SC                 : structural connectivity matrix (360 x 360)
%   RLFP_low_mean      : coupled intrinsic timescale (360 x 1)
%   RLFP_high_mean     : decoupled intrinsic timescale (360 x 1)
%   region_volume      : region volumes (360 x 1)
%   Nperm              : number of permutations (e.g., 1000)

%% ===============================
% 1. Observed correlations
% ===============================
SC_strength = sum(SC, 2);
region_volume = roi_volume;

[r_obs_low,  ~] = partialcorr(SC_strength, RLFP_low_mean,  region_volume, 'Type','Spearman');
[r_obs_high, ~] = partialcorr(SC_strength, RLFP_high_mean, region_volume, 'Type','Spearman');

%% ===============================
% 2. Null distributions
% ===============================
r_null_low  = zeros(Nperm,1);
r_null_high = zeros(Nperm,1);

upper_idx = find(triu(ones(size(SC)),1));
SC_vec = SC(upper_idx);

for p = 1:Nperm
    
    % ---- Randomize SC weights (preserve symmetry)
    SC_perm = zeros(size(SC));
    SC_perm(upper_idx) = SC_vec(randperm(length(SC_vec)));
    SC_perm = SC_perm + SC_perm';
    
    SC_strength_perm = sum(SC_perm,2);
    
    % ---- Partial correlations
    r_null_low(p) = partialcorr( ...
        SC_strength_perm, RLFP_low_mean, region_volume, 'Type','Spearman');
    
    r_null_high(p) = partialcorr( ...
        SC_strength_perm, RLFP_high_mean, region_volume, 'Type','Spearman');
end

%% ===============================
% 3. Z-scores and p-values
% ===============================
mu_low  = mean(r_null_low);   sigma_low  = std(r_null_low);
mu_high = mean(r_null_high);  sigma_high = std(r_null_high);

Z_low  = (r_obs_low  - mu_low)  / sigma_low;
Z_high = (r_obs_high - mu_high) / sigma_high;

p_low  = 2 * (1 - normcdf(abs(Z_low)));
p_high = 2 * (1 - normcdf(abs(Z_high)));

%% ===============================
% 4. Visualization (IEEE style)
% ===============================
figure('Color','w','Position',[150 200 1300 450])

% ---- Panel A: Coupled (low-pass)
subplot(1,2,1)
histogram(r_null_low,40,'Normalization','pdf', ...
    'FaceColor',[0.85 0.85 0.85],'EdgeColor','none'); hold on
xline(r_obs_low,'b','LineWidth',3)
xlabel('Partial Spearman r','FontSize',11)
ylabel('Probability density','FontSize',11)
title('Coupled (low-pass)','FontSize',12)
legend({'Null','Observed'},'Location','NorthWest')
box off

% ---- Panel B: Decoupled (high-pass)
subplot(1,2,2)
histogram(r_null_high,40,'Normalization','pdf', ...
    'FaceColor',[0.85 0.85 0.85],'EdgeColor','none'); hold on
xline(r_obs_high,'r','LineWidth',3)
xlabel('Partial Spearman r','FontSize',11)
title('Decoupled (high-pass)','FontSize',12)
legend({'Null','Observed'},'Location','NorthWest')
box off

sgtitle('Null model validation of structure–timescale coupling','FontSize',13)

%% ===============================
% 5. Console output (for Results)
% ===============================
fprintf('\nFIGURE 6 — Null model validation\n');

fprintf('Coupled (low-pass):\n');
fprintf('  Observed r = %.4f\n', r_obs_low);
fprintf('  Null mean = %.4f ± %.4f\n', mu_low, sigma_low);
fprintf('  Z = %.2f, p = %.4g\n\n', Z_low, p_low);

fprintf('Decoupled (high-pass):\n');
fprintf('  Observed r = %.4f\n', r_obs_high);
fprintf('  Null mean = %.4f ± %.4f\n', mu_high, sigma_high);
fprintf('  Z = %.2f, p = %.4g\n', Z_high, p_high);

end
