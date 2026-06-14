function plot_Figure3_Replication_H2()
% ============================================================
% Figure 3
% Replication of structure–timescale relationship (H2)
% Based on Fulcher et al., using HCP data
%
% Panels:
%   (A) Scatter: SC strength vs RLFP (original BOLD)
%   (B) Partial Spearman correlation (group-level)
%   (C) Spatial cortical map of RLFP
%
% Assumes variables exist in base workspace:
%   SC_strength
%   RLFP_mean
%   RegionVolume
%   Glasser_labels (optional)
% ============================================================

%% ------------------------------------------------------------
% Load data from base workspace
% ------------------------------------------------------------
SC_strength  = evalin('base','SC_strength');
RLFP_mean    = evalin('base','RLFP_mean');
RegionVolume = evalin('base','roi_volume');

%% ------------------------------------------------------------
% Partial Spearman correlation (control for volume)
% ------------------------------------------------------------
valid_idx = ...
    ~isnan(SC_strength) & ...
    ~isnan(RegionVolume);
[r_H2, p_H2] = partialcorr( ...
    SC_strength(valid_idx), RLFP_mean(valid_idx), RegionVolume(valid_idx), 'Type','Spearman');

%% ------------------------------------------------------------
% Create figure
% ------------------------------------------------------------
figure('Color','w','Position',[150 150 1400 450])

%% ============================================================
% (A) Scatter: SC vs RLFP
% ============================================================
subplot(1,3,1)
scatter(SC_strength(valid_idx), RLFP_mean(valid_idx), 35, 'filled')
lsline
xlabel('Structural connectivity strength')
ylabel('RLFP (original BOLD)')
title('Region-wise structure–timescale relationship')
grid on

text(min(SC_strength(valid_idx)), max(RLFP_mean(valid_idx))*0.95, ...
    sprintf('r = %.2f, p = %.3g', r_H2, p_H2), ...
    'FontSize',10,'FontWeight','bold');

%% ============================================================
% (B) Partial correlation summary
% ============================================================
subplot(1,3,2)
bar(r_H2)
set(gca,'XTick',1,'XTickLabel',{'All regions'},'FontSize',11)
ylabel('Partial Spearman r')
title('Group-level coupling')
yline(0,'--')
grid on

text(1, r_H2*0.9, sprintf('p = %.3g', p_H2), ...
    'HorizontalAlignment','center','FontSize',10);

end
