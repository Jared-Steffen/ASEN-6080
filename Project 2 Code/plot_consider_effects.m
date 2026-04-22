function plot_consider_effects(t,P,~,S,Pcc,consider_labels)

[~,~,N] = size(P);
consider_labels = consider_labels(:).';

% ---- Neon colors ----
neonColors = [ ...
    0.0  1.0  1.0;   % cyan
    1.0  0.0  1.0;   % magenta
    0.0  1.0  0.3;   % green
    1.0  0.0  0.3;   % red
    1.0  0.3  0.0;   % orange
    0.3  0.6  1.0;   % blue
    1.0  1.0  0.0;   % yellow
    0.6  0.0  1.0;   % purple
    0.0  1.0  0.6;   % mint
    0.0  0.8  1.0;   % sky
    0.6  1.0  0.0;   % lime
    1.0  0.6  0.0];  % amber

% Groups of consider parameters
group_idx = { ...
    [1 2], ...          % muE, muS
    [3 4 5], ...        % station 34
    [6 7 8], ...        % station 65
    [9 10 11]};         % station 13

group_names = { ...
    'Gravitational Parmeters', ...
    'DSS 34 Position', ...
    'DSS 65 Position', ...
    'DSS 13 Position'};

state_groups = {1:3, 4:6, 7};
state_group_names = {'Position','Velocity','C_R'};
state_labels_pos = {'x','y','z'};
state_labels_vel = {'v_x','v_y','v_z'};

for g = 1:length(group_idx)

    idx_consider = group_idx{g};

    %% -------- POSITION --------
    figure;
    for i = 1:3
        subplot(3,1,i); hold on;

        h = [];
        legend_entries = {};

        % Nominal
        sig_nom = zeros(N,1);
        for k = 1:N
            sig_nom(k) = 3*sqrt(P(i,i,k));
        end

        c = neonColors(1,:);
        line = plot(t/3600, +sig_nom, '-', 'Color', c);
        plot(t/3600, -sig_nom, '-', 'Color', c, 'HandleVisibility','off');

        h = [h, line];
        legend_entries{end+1} = 'Nominal';

        % Selected consider params only
        for m = 1:length(idx_consider)
            j = idx_consider(m);

            sig_j = zeros(N,1);
            for k = 1:N
                sj = S(:,j,k);
                Pxx_j = P(:,:,k) + sj*Pcc(j,j)*sj.';
                sig_j(k) = 3*sqrt(Pxx_j(i,i));
            end

            c = neonColors(m+1,:);
            line = plot(t/3600, +sig_j, '--', 'Color', c);
            plot(t/3600, -sig_j, '--', 'Color', c, 'HandleVisibility','off');

            h = [h, line];
            legend_entries{end+1} = char(consider_labels(j));
        end

        title([group_names{g} ' - ' state_labels_pos{i} ' \pm 3\sigma'])
        xlabel('Time [hr]')
        ylabel('\pm 3\sigma [km]')
        grid on
        lgd = legend(h, legend_entries, 'Location','best');
        lgd.Position = [0.75 0.935 0.25 0.05];
    end

    %% -------- VELOCITY --------
    figure;
    for i = 4:6
        subplot(3,1,i-3); hold on;

        h = [];
        legend_entries = {};

        % Nominal
        sig_nom = zeros(N,1);
        for k = 1:N
            sig_nom(k) = 3*sqrt(P(i,i,k));
        end

        c = neonColors(1,:);
        line = plot(t/3600, +sig_nom, '-', 'Color', c);
        plot(t/3600, -sig_nom, '-', 'Color', c, 'HandleVisibility','off');

        h = [h, line];
        legend_entries{end+1} = 'Nominal';

        % Selected consider params only
        for m = 1:length(idx_consider)
            j = idx_consider(m);

            sig_j = zeros(N,1);
            for k = 1:N
                sj = S(:,j,k);
                Pxx_j = P(:,:,k) + sj*Pcc(j,j)*sj.';
                sig_j(k) = 3*sqrt(Pxx_j(i,i));
            end

            c = neonColors(m+1,:);
            line = plot(t/3600, +sig_j, '--', 'Color', c);
            plot(t/3600, -sig_j, '--', 'Color', c, 'HandleVisibility','off');

            h = [h, line];
            legend_entries{end+1} = char(consider_labels(j));
        end

        title([group_names{g} ' - ' state_labels_vel{i-3} ' \pm 3\sigma'])
        xlabel('Time [hr]')
        ylabel('\pm 3\sigma [km/s]')
        grid on
        lgd = legend(h, legend_entries, 'Location','best');
        lgd.Position = [0.75 0.935 0.25 0.05];
    end

    %% -------- CR --------
    figure; hold on;

    i = 7;
    h = [];
    legend_entries = {};

    % Nominal
    sig_nom = zeros(N,1);
    for k = 1:N
        sig_nom(k) = 3*sqrt(P(i,i,k));
    end

    c = neonColors(1,:);
    line = plot(t/3600, +sig_nom, '-', 'Color', c);
    plot(t/3600, -sig_nom, '-', 'Color', c, 'HandleVisibility','off');

    h = [h, line];
    legend_entries{end+1} = 'Nominal';

    % Selected consider params only
    for m = 1:length(idx_consider)
        j = idx_consider(m);

        sig_j = zeros(N,1);
        for k = 1:N
            sj = S(:,j,k);
            Pxx_j = P(:,:,k) + sj*Pcc(j,j)*sj.';
            sig_j(k) = 3*sqrt(Pxx_j(i,i));
        end

        c = neonColors(m+1,:);
        line = plot(t/3600, +sig_j, '--', 'Color', c);
        plot(t/3600, -sig_j, '--', 'Color', c, 'HandleVisibility','off');

        h = [h, line];
        legend_entries{end+1} = char(consider_labels(j));
    end

    title([group_names{g} ' - C_R \pm 3\sigma'])
    xlabel('Time [hr]')
    ylabel('\pm 3\sigma')
    grid on
    lgd = legend(h, legend_entries, 'Location','best');
    lgd.Position = [0.75 0.935 0.25 0.05];

end

end
