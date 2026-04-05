function twoSigmaCrossStats(Xmc, Xnom, P, t, time_idx, labels)
% TWOSIGMACROSSSTATS
% Reports percentage of Monte Carlo runs outside:
%   1) componentwise 2-sigma bounds
%   2) pairwise 2-sigma covariance ellipses
%
% Inputs:
%   Xmc      : Nt x d x Nmc Monte Carlo state history
%   Xnom     : Nt x d nominal / propagated estimate
%   P        : d x d x Nt propagated covariance history
%   t        : time vector
%   time_idx : time index to evaluate
%   labels   : optional cell array of state labels

    if nargin < 6 || isempty(labels)
        d = size(Xmc,2);
        labels = arrayfun(@(i) sprintf('x_%d',i), 1:d, 'UniformOutput', false);
    end

    % MC samples at this time: Nmc x d
    X = squeeze(Xmc(time_idx,:,:))';

    % Nominal and covariance at this time
    x_nom = Xnom(time_idx,:);
    Pk = P(:,:,time_idx);

    % ------------------------------------------------------------
    % 1D componentwise results
    % ------------------------------------------------------------
    sigma = sqrt(diag(Pk))';
    outside_each = abs(X - x_nom) > 2*sigma;
    pct_outside_each = 100 * sum(outside_each,1) / size(X,1);

    fprintf('\n============================================================\n');
    fprintf('2-sigma outside percentages at t = %.2f hours\n', t(time_idx)/3600);
    fprintf('============================================================\n');

    fprintf('\nComponentwise 1D bounds:\n');
    for i = 1:length(labels)
        fprintf('%-3s | Outside = %7.2f%%\n', labels{i}, pct_outside_each(i));
    end

    % ------------------------------------------------------------
    % 2D pairwise ellipse results
    % ------------------------------------------------------------
    fprintf('\nPairwise 2D covariance ellipses:\n');

    dstate = size(X,2);
    pct_outside_pairs = nan(dstate,dstate);

    for i = 2:dstate
        for j = 1:i-1
            X2 = X(:,[j i]);                 % Nmc x 2
            mu2 = [x_nom(j); x_nom(i)];      % 2 x 1
            P2 = Pk([j i],[j i]);            % 2 x 2

            outside_pair = false(size(X2,1),1);

            for k = 1:size(X2,1)
                dx = X2(k,:)' - mu2;
                d2 = dx' * (P2 \ dx);
                outside_pair(k) = d2 > 4;    % outside 2-sigma ellipse
            end

            pct_outside_pairs(i,j) = 100 * sum(outside_pair) / size(X2,1);

            fprintf('%-3s vs %-3s | Outside = %7.2f%%\n', ...
                labels{j}, labels{i}, pct_outside_pairs(i,j));
        end
    end
end