function plot_Figure5_InterIndividual_H4()
% ============================================================
% Figure 5
% Inter-individual robustness of structure–timescale coupling
% Tests Hypothesis H4
%
% Panels:
% (A) Subject-wise partial correlations (coupled vs decoupled)
% (B) Paired boxplot comparison
% (C) Effect size (Cohen's d, paired)
%
% Required variables in base workspace:
%   SC_strength_subj   [360 x Nsub]
%   RLFP_low_subj     [360 x Nsub]
%   RLFP_high_subj    [360 x Nsub]
%   RegionVolume      [360 x 1]
% ============================================================

%% ------------------------------------------------------------
% Load variables
% ------------------------------------------------------------
SC_all   = evalin('base','SC_strength');
RLFP_low_all  = evalin('base','RLFP_low');
RLFP_high_all = evalin('base','RLFP_high');
RegionVolume  = evalin('base','roi_volume');

% Nsub = size(SC_all,2);
Nsub = size(RLFP_low_all,2);

r_low_subj  = zeros(Nsub,1);
r_high_subj = zeros(Nsub,1);

%% ------------------------------------------------------------
% Subject-wise partial correlations
% ------------------------------------------------------------
for s = 1:Nsub
    
    SC_s   = SC_all;
    low_s  = RLFP_low_all(:,s);
    high_s = RLFP_high_all(:,s);
    
    r_low_subj(s) = partialcorr(SC_s, low_s,  RegionVolume, 'Type','Spearman');
    r_high_subj(s)= partialcorr(SC_s, high_s, RegionVolume, 'Type','Spearman');
end

%% ------------------------------------------------------------
% Paired statistical test
% ------------------------------------------------------------
[p_wilcoxon,~,stats] = signrank(r_low_subj, r_high_subj);

%% ------------------------------------------------------------
% Effect size: Cohen's d (paired)
% ------------------------------------------------------------
diff_r = r_low_subj - r_high_subj;
cohen_d = mean(diff_r) / std(diff_r);

%% ------------------------------------------------------------
% Create figure
% ------------------------------------------------------------
figure('Color','w','Position',[150 150 1500 450])

%% ============================================================
% (A) Subject-wise r values
% ============================================================
subplot(1,3,1)
plot(r_low_subj, r_high_subj, 'o','MarkerFaceColor',[0.2 0.4 0.8])
hold on
plot([-0.5 0.5],[-0.5 0.5],'k--')
xlabel('r (Coupled)')
ylabel('r (Decoupled)')
title('Subject-wise structure–timescale coupling')
axis square
grid on

%% ============================================================
% (B) Paired boxplot
% ============================================================
subplot(1,3,2)
boxplot([r_low_subj, r_high_subj], ...
        'Labels',{'Coupled','Decoupled'})
ylabel('Partial Spearman r')
title('Paired comparison across subjects')
grid on

text(1.5, max([r_low_subj; r_high_subj])*0.9, ...
    sprintf('Wilcoxon p = %.3g', p_wilcoxon), ...
    'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');

%% ============================================================
% (C) Effect size
% ============================================================
subplot(1,3,3)
bar(cohen_d)
set(gca,'XTick',1,'XTickLabel',{'Cohen''s d'})
ylabel('Effect size')
title('Magnitude of coupled > decoupled effect')
yline(0,'--')
grid on

text(1, cohen_d*0.9, sprintf('d = %.2f', cohen_d), ...
    'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');

sgtitle('Figure 5 | Inter-individual robustness of structure–timescale coupling', ...
        'FontWeight','bold');

end
