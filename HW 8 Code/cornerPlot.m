function cornerPlot(Xmc,Xnom,t,time_idx,labels)
% CORNER_PLOT_MC  Corner plot from Monte Carlo state histories
%
% Inputs:
%   Xmc      : Monte Carlo results
%              size = Nt x d x Nmc
%   time_idx : time index to plot
%   labels   : optional cell array of state labels
%
% Example:
%   corner_plot_mc(Xstate, 1, {'x','y','z','vx','vy','vz'}, 25)
%
%   % plot at 12 hr if t is your time vector:
%   [~,k] = min(abs(t - 12*3600));
%   corner_plot_mc(Xstate, k, {'x','y','z','vx','vy','vz'}, 30)

    if nargin < 3 || isempty(labels)
        d = size(Xmc,2);
        labels = arrayfun(@(i) sprintf('x_%d',i), 1:d, 'UniformOutput', false);
    end

    [Nt, d, Nmc] = size(Xmc);


    % Extract samples at requested time:
    % X = Nmc x d
    X = squeeze(Xmc(time_idx,:,:))';

    % Nominal at this time
    x_nom = Xnom(time_idx,:);

    % Statistics
    mu = mean(X,1);
    sigma = std(X,0,1);

    % ============================================================
    % Print
    % ============================================================
    fprintf('\n========================================\n');
    fprintf('Statistics at t = %.2f hours\n', t(time_idx)/3600);
    fprintf('========================================\n');

    for i = 1:d
        fprintf('%-3s | Nominal = %12.6f | Mean = %12.6f | Sigma = %12.6f\n', ...
                labels{i}, x_nom(i), mu(i), sigma(i));
    end

    % ============================================================
    % Plot
    % ============================================================
    figure;
    for i = 1:d
        for j = 1:d
            subplot(d,d,(i-1)*d + j)

            if i == j
                histogram(X(:,i), 'Normalization', 'pdf');
                hold on
                xline(mu(i), 'r', 'LineWidth', 1.2);
                xline(x_nom(i), 'k--', 'LineWidth', 1.2);
                grid on

            elseif i > j
                mu2 = [mu(j); mu(i)];
                X2 = X(:,[j i]);
                P2 = cov(X2);
            
                [V,D] = eig(P2);
                th = linspace(0, 2*pi, 200);
                circle = [cos(th); sin(th)];
            
                hold on
            
                % --- Monte Carlo samples ---
                scatter(X2(:,1), X2(:,2), 8, 'filled', ...
                    'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
            
                % --- Get default color order ---
                co = get(gca,'ColorOrder');
                scales = [1 2 3];
                styles = {'-','--',':'};
            
                for s = 1:3
                    ellipse = scales(s) * V * sqrt(D) * circle;
                    pts = ellipse + mu2;
            
                    plot(pts(1,:), pts(2,:), ...
                        'LineStyle', styles{s}, ...
                        'LineWidth', 1.5, ...
                        'Color', co(s,:));
                end
            
                % --- Mean ---
                plot(mu(j), mu(i), 'o', 'LineWidth', 1.5);
            
                % --- Nominal ---
                plot(x_nom(j), x_nom(i), 'x', 'LineWidth', 1.8);
            
                grid on
            else
                axis off
                continue
            end

            if i == d
                xlabel(labels{j})
            else
                set(gca,'XTickLabel',[])
            end

            if j == 1
                ylabel(labels{i})
            else
                set(gca,'YTickLabel',[])
            end
        end
    end
end