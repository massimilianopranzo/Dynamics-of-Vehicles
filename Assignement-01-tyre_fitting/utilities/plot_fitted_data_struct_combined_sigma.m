function [] = plot_fitted_data_struct_combined_sigma(x_raw, y_raw, x_fit, y_fit, x_fit_sig, y_fit_sig, ...
  label_x, label_y, data_label, leg_angle,legend_names, name, plot_title, line_width, ...
  font_size_title, colors_vector)

	if ismatrix(data_label) == 0;
		error('Data label must be a vector')
		return;
	end


	len = length(data_label);
	rows = ceil(len/2);
	
	fig = figure('Color','w');
  for i = 1:len
% 		if i == len && mod(len, 2) == 1
% 			subplot(rows, 2, i + 0.5)
% 		else
			subplot(rows, 2, i)
% 		end
		grid on
		hold on

		markers = ['.', 'o', '+'];
		for j = 1:length(leg_angle) %loop over alpha
      		plot(x_fit(:), y_fit{i}{j}(:),'-', 'DisplayName', leg_angle(j) , 'LineWidth', line_width, 'Color', colors_vector(j,:))
      		plot(x_fit_sig(:), y_fit_sig{i}{j}(:),'--', 'DisplayName', leg_angle(j) , 'LineWidth', line_width, 'Color', colors_vector(j,:))
			if ~isempty(x_raw) && ~isempty(y_raw)
			    plot(x_raw{i}{j}(:), y_raw{i}{j}(:), markers(j), 'DisplayName','Raw', 'LineWidth', .2 * line_width, 'Color', colors_vector(j,:)) %  colors_vector(j,:)
			end

		end
		% [x_fit{i}, order] = sort(x_fit{i});
		% y_fit{i} = y_fit{i}(order);
		title([data_label(i)], 'interpreter','latex', 'FontSize', font_size_title);
		xlabel(label_x)
		ylabel(label_y)
		hold off
  end
  if ~isempty(legend_names)
  legend(legend_names, 'position', [0.679258117202733,0.238193265647898,0.127430735291504,0.058136925199618], 'FontSize', font_size_title,'NumColumns',2 )
  end
	sgtitle(plot_title, 'interpreter','latex', 'FontSize', font_size_title)
	
	fig_name = ['images\', name, '.png'];
	export_fig(fig, fig_name)
	
end