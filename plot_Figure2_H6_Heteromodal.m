function plot_Figure2_H6_Heteromodal()
% ============================================================
% Figure 2
% Structure–timescale coupling of decoupled dynamics
% Heteromodal vs Non-heteromodal cortex
%
% Assumes variables exist in base workspace:
%   SC_strength
%   RLFP_high_mean
%   RegionVolume
%   idx_heteromodal
% ============================================================

%% ------------------------------------------------------------
% Load variables from base workspace
% ------------------------------------------------------------
SC_strength    = evalin('base','SC_strength');
RLFP_high_mean = evalin('base','RLFP_high_mean');
RegionVolume   = evalin('base','roi_volume');
idx_hetero     = evalin('base','idx_heteromodal');
Z_diff         = evalin('base','z_diff'); 
p_diff         = evalin('base','p_diff'); 

% Define non-heteromodal regions
idx_nonhetero = setdiff(1:length(SC_strength), idx_hetero);

%% ------------------------------------------------------------
% Partial correlations (controlling for region volume)
% ------------------------------------------------------------
% [r_hetero, ~] = partialcorr( ...
%     SC_strength(idx_hetero), ...
%     RLFP_high_mean(idx_hetero), ...
%     RegionVolume(idx_hetero), ...
%     'Type','Spearman');
% 
% [r_nonhetero, ~] = partialcorr( ...
%     SC_strength(idx_nonhetero), ...
%     RLFP_high_mean(idx_nonhetero), ...
%     RegionVolume(idx_nonhetero), ...
%     'Type','Spearman');
[r_hetero, p_hetero] = partialcorr( ...
    SC_strength(idx_hetero), ...
    RLFP_high_mean(idx_hetero), ...
    RegionVolume(idx_hetero), ...
    'Type','Spearman');

[r_nonhetero, p_nonhetero] = partialcorr( ...
    SC_strength(idx_nonhetero), ...
    RLFP_high_mean(idx_nonhetero), ...
    RegionVolume(idx_nonhetero), ...
    'Type','Spearman');

%% ------------------------------------------------------------
% Mean RLFP_high
% ------------------------------------------------------------
mean_RLFP = [ ...
    mean(RLFP_high_mean(idx_hetero)), ...
    mean(RLFP_high_mean(idx_nonhetero)) ];

%% ============================================================
% Plot Figure 2
% ============================================================
figure('Color','w','Position',[200 200 1400 450])


%% ---------------------------
% (A) Partial correlation bar
% ---------------------------
subplot(1,4,1)
bar([r_hetero, r_nonhetero])
set(gca,'XTickLabel',{'Transmodal','Unimodal'})
ylabel('Partial Spearman r')
yline(0,'--')
grid on

text(1, r_hetero, sprintf('r=%.2f\np=%.3f', r_hetero, p_hetero), ...
     'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);

text(2, r_nonhetero, sprintf('r=%.2f\np=%.3f', r_nonhetero, p_nonhetero), ...
     'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);
%% ---------------------------
% (B) Mean RLFP_high bar
% ---------------------------
subplot(1,4,2)
bar(mean_RLFP)
set(gca,'XTickLabel',{'Transmodal','Unimodal'}, ...
        'FontSize',11)
ylabel('Mean RLFP_{high}')
title('Decoupled timescales')
grid on
%% ------------------------------------------------------------
% Annotate Fisher Z-test result (subplot 1,4,2)
% ------------------------------------------------------------

subplot(1,4,2)   % ensure correct panel is active
ylim_current = ylim;   % get y-axis limits

text(1.5, ylim_current(2)*0.85, ...
    sprintf('Fisher Z = %.2f\np = %.4f', Z_diff, p_diff), ...
    'HorizontalAlignment','center', ...
    'FontSize',10, ...
    'FontWeight','bold');

%% ---------------------------
% (C) Scatter: heteromodal
% ---------------------------
subplot(1,4,3)
hold on
scatter(SC_strength(idx_hetero), ...
        RLFP_high_mean(idx_hetero), ...
        35,'filled')
lsline
xlabel('SC strength')
ylabel('RLFP_{high}')
title('Transmodal cortex')
grid on
hold off

sgtitle('Figure 2 | Structure–timescale coupling of decoupled dynamics', ...
        'FontWeight','bold')
%% ---------------------------
% (D) Scatter: non-heteromodal
% ---------------------------
subplot(1,4,4)
hold on
scatter(SC_strength(idx_nonhetero), ...
        RLFP_high_mean(idx_nonhetero), ...
        35,'filled')
lsline
xlabel('SC strength')
ylabel('RLFP_{high}')
title('Unimodal cortex')
grid on
hold off


end
