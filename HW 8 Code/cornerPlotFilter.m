function cornerPlotFilter(Xmc,Xnom,P,t,time_idx,labels)
% CORNERPLOT  Corner plot using Monte Carlo dots and propagated covariance ellipses
%
% Inputs:
%   Xmc      : Monte Carlo results, size = Nt x d x Nmc
%   Xnom     : nominal / propagated estimate, size = Nt x d
%   P        : propagated covariance history, size = d x d x Nt
%   t        : time vector
%   time_idx : time index to plot
%   labels   : optional cell array of state labels

    if nargin < 6 || isempty(labels)
        d = size(Xmc,2);
        labels = arrayfun(@(i) sprintf('x_%d',i), 1:d, 'UniformOutput', false);
    end

    [Nt, d, ~] = size(Xmc);

    % Extract MC samples at requested time: X = Nmc x d
    X = squeeze(Xmc(time_idx,:,:))';

    % Nominal / propagated estimate at this time
    x_nom = Xnom(time_idx,:);

    % Sample mean from MC
    mu_mc = mean(X,1);

    % Sigma from propagated covariance
    sigma_ckf = sqrt(diag(P(:,:,time_idx)))';

    % ============================================================
    % Print
    % ============================================================
    fprintf('\n========================================\n');
    fprintf('Statistics at t = %.2f hours\n', t(time_idx)/3600);
    fprintf('========================================\n');

    for i = 1:d
        fprintf('%-3s | Nominal = %12.6f | Mean = %12.6f | Sigma = %12.6f\n', ...
                labels{i}, x_nom(i), mu_mc(i), sigma_ckf(i));
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
                xline(mu_mc(i), 'LineWidth', 1.2);
                xline(x_nom(i), '--', 'LineWidth', 1.2);

                % 1σ, 2σ, 3σ bounds from propagated covariance
                co = get(gca,'ColorOrder');
                scales = [1 2 3];
                styles = {'-','--',':'};

                for s = 1:3
                    xline(x_nom(i) + scales(s)*sigma_ckf(i), ...
                        'LineStyle', styles{s}, ...
                        'LineWidth', 1.2, ...
                        'Color', co(s,:));
                    xline(x_nom(i) - scales(s)*sigma_ckf(i), ...
                        'LineStyle', styles{s}, ...
                        'LineWidth', 1.2, ...
                        'Color', co(s,:));
                end

                grid on

            elseif i > j
                X2 = X(:,[j i]);
                mu2 = [x_nom(j); x_nom(i)]; 
                P2 = P([j i],[j i],time_idx);   

                [V,D] = eig(P2);
                th = linspace(0, 2*pi, 200);
                circle = [cos(th); sin(th)];

                hold on

                % Monte Carlo samples as dots
                scatter(X2(:,1), X2(:,2), 8, 'filled', ...
                    'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);

                % Default color order
                co = get(gca,'ColorOrder');
                scales = [1 2 3];
                styles = {'-','--',':'};

                % 1σ, 2σ, 3σ ellipses from propagated covariance
                for s = 1:3
                    ellipse = scales(s) * V * sqrt(D) * circle;
                    pts = ellipse + mu2;

                    plot(pts(1,:), pts(2,:), ...
                        'LineStyle', styles{s}, ...
                        'LineWidth', 1.5, ...
                        'Color', co(s,:));
                end

                % MC mean
                plot(mu_mc(j), mu_mc(i), 'o', 'LineWidth', 1.5);

                % Nominal / propagated estimate
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