function [] = plot_fitted_data(x_raw, y_raw, x_fit, y_fit, label_x, label_y, data_label, name, plot_title, line_width, font_size_title)

	if ismatrix(data_label) == 0
		error('Data label must be a vector')
		return;
	end
  fig = figure('Color','w');
	plot(x_raw, y_raw,'.', 'DisplayName','Raw data', 'LineWidth', line_width)
	hold on
	for i = 1:length(data_label)
			plot(x_fit(:,i), y_fit(:,i),'-', 'DisplayName', data_label(i), 'LineWidth', line_width)
	end
	grid on
	legend('Location', 'best')
	xlabel(label_x)
	ylabel(label_y)
	title(plot_title, 'interpreter','latex', 'FontSize', font_size_title)
	
	fig_name = ['images\', name, '.svg'];
	export_fig(fig, fig_name)
	
end