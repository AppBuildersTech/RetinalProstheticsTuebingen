function T02_plots(STA_ps, D_ps, exp_ps)

line_thickness = 2;

estim_meanline = STA_ps.estim_mean * ones(2 * exp_ps.tKerLen, 1);

if exp_ps.Normalize == 1
    plt_ylim = [-1, 1];
else
	plt_ylim = [-1200, -400];
end
yaxis_line = zeros(length(plt_ylim(1):100:plt_ylim(2)));

fig_basename = sprintf('%s_[%s]',exp_ps.exp_id,exp_ps.cell_id);

%% Plot 1: STA_ps.STA
figure
set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
set(gcf, 'color', 'w');

plot(STA_ps.STA_t, STA_ps.STA,'LineWidth', line_thickness);hold on;
plot(STA_ps.STA_t, estim_meanline, 'k');

plot(yaxis_line, plt_ylim(1):100:plt_ylim(2), 'k');

title(fig_basename, 'Interpreter', 'none')
ylim([plt_ylim(1) plt_ylim(2)])

%saveas(gcf, [exp_ps.work_dir, fig_basename,'_STA.fig']);
saveas(gcf, [exp_ps.work_dir, fig_basename, '_STA.jpeg']);

end