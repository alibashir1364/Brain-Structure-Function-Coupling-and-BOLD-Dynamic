function plot_Figure4_Coupled_vs_Decoupled()
% ============================================================
% Figure 4
% Coupled vs Decoupled intrinsic timescales (Group-level)
% Tests Hypotheses H1 and H3
%
% Panels:
% (A) SC vs RLFP_low  (coupled / low-pass)
% (B) SC vs RLFP_high (decoupled / high-pass)
% (C) Partial correlation comparison
% (D) Cortical difference map (low - high)
%
% Required variables in base workspace:
%   SC_strength
%   RLFP_low_mean
%   RLFP_high_mean
%   RegionVolume
% ============================================================

%% ------------------------------------------------------------
% Load variables
% ------------------------------------------------------------
SC_strength    = evalin('base','SC_strength');
RLFP_low       = evalin('base','RLFP_low_mean');
RLFP_high      = evalin('base','RLFP_high_mean');
RegionVolume   = evalin('base','roi_volume');

%% ------------------------------------------------------------
% Partial correlations (Spearman, volume-controlled)
% ------------------------------------------------------------
[r_low,  p_low]  = partialcorr(SC_strength, RLFP_low,  RegionVolume, 'Type','Spearman');
[r_high, p_high] = partialcorr(SC_strength, RLFP_high, RegionVolume, 'Type','Spearman');

%% ------------------------------------------------------------
% Fisher Z-test: low vs high
% ------------------------------------------------------------
Z_low  = atanh(r_low);
Z_high = atanh(r_high);

N = length(SC_strength);
Z_diff = (Z_low - Z_high) / sqrt(2/(N-3));
p_diff = 2 * (1 - normcdf(abs(Z_diff)));

%% ------------------------------------------------------------
% Create figure
% ------------------------------------------------------------
figure('Color','w','Position',[120 120 1500 450])

%% ============================================================
% (A) Scatter: SC vs RLFP_low
% ============================================================
subplot(1,4,1)
scatter(SC_strength, RLFP_low, 30, 'filled')
lsline
xlabel('Structural connectivity strength')
ylabel('RLFP_{low} (coupled)')
title('Coupled dynamics')
grid on

text(min(SC_strength), max(RLFP_low)*0.95, ...
    sprintf('r = %.2f\np = %.3g', r_low, p_low), ...
    'FontSize',10,'FontWeight','bold');

%% ============================================================
% (B) Scatter: SC vs RLFP_high
% ============================================================
subplot(1,4,2)
scatter(SC_strength, RLFP_high, 30, 'filled')
lsline
xlabel('Structural connectivity strength')
ylabel('RLFP_{high} (decoupled)')
title('Decoupled dynamics')
grid on

text(min(SC_strength), max(RLFP_high)*0.95, ...
    sprintf('r = %.2f\np = %.3g', r_high, p_high), ...
    'FontSize',10,'FontWeight','bold');

%% ============================================================
% (C) Partial correlation comparison
% ============================================================
subplot(1,4,3)
bar([r_low, r_high])
set(gca,'XTickLabel',{'Coupled','Decoupled'},'FontSize',11)
ylabel('Partial Spearman r')
title('Structure–timescale coupling')
yline(0,'--')
grid on

ylim_current = ylim;
text(1.5, ylim_current(2)*0.85, ...
    sprintf('Fisher Z = %.2f\np = %.3g', Z_diff, p_diff), ...
    'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');

%% ============================================================
% (D) Cortical difference map (low - high)
% ============================================================
subplot(1,4,4)

diff_map = RLFP_low_mean - RLFP_high_mean;

if evalin('base','exist(''plot_glasser_surface'',''file'')')
    evalin('base','plot_glasser_surface(diff_map)')
    title('\Delta Timescale (low - high)')
else
    imagesc(diff_map')
    colormap(parula)
    colorbar
    title('\Delta RLFP (low - high)')
    xlabel('Region index (Glasser)')
end

sgtitle('Figure 4 | Structural connectivity selectively constrains graph-smooth activity', ...
        'FontWeight','bold');

end
