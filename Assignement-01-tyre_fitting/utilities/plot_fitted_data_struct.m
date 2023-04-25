function [] = plot_fitted_data_struct(x_raw, y_raw, x_fit, y_fit, ...
  label_x, label_y, data_label, name, plot_title, line_width, ...
  font_size_title, colors_vector)

	if ismatrix(data_label) == 0;
		error('Data label must be a vector')
		return;
	end


	len = length(data_label);
	rows = ceil(len/2);
	
	fig = figure('Color','w');
  for i = 1:len
		if i == len && mod(len, 2) == 1
			subplot(rows, 2, i + 0.5)
		else
			subplot(rows, 2, i)
		end
		grid on
		hold on
		plot(x_raw{i}, y_raw{i},'.', 'DisplayName','Raw', 'LineWidth', line_width, 'Color', colors_vector(1,:))
		[x_fit{i}, order] = sort(x_fit{i});
		y_fit{i} = y_fit{i}(order);
		plot(x_fit{i}, y_fit{i},'-', 'DisplayName', 'Fit' , 'LineWidth', line_width, 'Color', colors_vector(2,:))
		title([data_label(i)], 'interpreter','latex', 'FontSize', font_size_title);
    xlabel(label_x)
    ylabel(label_y)
		hold off
  end
	legend('Location', 'best')
	sgtitle(plot_title, 'interpreter','latex', 'FontSize', font_size_title)
	
	fig_name = ['images\', name, '.png'];
	export_fig(fig, fig_name)
	
end