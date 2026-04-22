function plot_Bplane(cov_BPlane_list,BdotRhat_list,BdotThat_list,DCO_t_list)

figure(); hold on; grid on;
xlabel('B \cdot T'); ylabel('B \cdot R');
title('B-Plane Analysis');

nCases = length(BdotRhat_list);
colors = lines(nCases);

legend_entries = cell(nCases,1);
h = gobjects(nCases,1);

for i = 1:nCases
    
    % Get TR covariance
    P = cov_BPlane_list(2:3,2:3,i);
    
    % Eigen decomposition
    [V,D] = eig(P);
    
    % 3-sigma ellipse
    t = linspace(0,2*pi,200);
    ellipse = 3 * V * sqrt(D) * [cos(t); sin(t)];
    
    % Center
    x0 = BdotThat_list(i);
    y0 = BdotRhat_list(i);
    
    % Plot ellipse
    h(i) = plot(x0 + ellipse(1,:), y0 + ellipse(2,:), ...
        'Color', colors(i,:), 'LineWidth', 1.5);
    
    % Plot point
    plot(x0, y0, '.', 'Color', colors(i,:), 'MarkerSize', 20);
    
    % Legend entry
    if iscell(DCO_t_list)
        val = DCO_t_list{i};
    else
        val = DCO_t_list(i);
    end

    if isnumeric(val)
        legend_entries{i} = sprintf('DCO = %.2f days', val);
    else
        legend_entries{i} = char(val);
    end
end

legend(h, legend_entries, 'Location', 'best');

end