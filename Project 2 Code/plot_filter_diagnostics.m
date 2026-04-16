function plot_filter_diagnostics(t,Xhat,P,y,yhat,R,station_id,filter_type,data_type)
%{
Inputs:
    >t: Nx1 time vector [s]
    >Xtrue: Nx6 true state history [km, km/s]
    >Xhat: Nx6 estimated state history [km, km/s]
    >P: 6x6xN state covariance history
    >y: Nx2 pre-fit residuals  [range(km), range-rate(km/s)] (NaNs allowed)
    >yhat: Nx2 post-fit residuals [range(km), range-rate(km/s)] (NaNs allowed)
    >R: measurement noise covariance matrix
    >station_id: Nx1 station ID associated with each residual time step
    >filter_type: string of filter name for plots
    >data_type: string of "range" or "range rate" if only considering
                one (or empty for both)
Outputs:
    > State errors vs time with ±3σ covariance bounds
    > Pre-fit and post-fit residuals (time series + QQ), color-coded by station
%}

%% Plot formatting
set(groot,'defaultFigureColor','w')
set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesLineWidth',1.2)
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesGridAlpha',0.3)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
lineStyles   = {'--','-',':','-.'};
markerStyles = {'.','o','x'};

% 3-sigma bounds
n = size(Xhat,2);
N = length(t);
sig = zeros(N,n);
for i = 1:N
    sig(i,:) = sqrt(diag(P(1:n,1:n,i))).';
    if ~isreal(sig(i,:))
        disp("Non-real sigma at index i = " + i)
    end
end
sig3 = 3*sig;


%% State error plots
r_labels = {'x [km]','y [km]','z [km]'};
v_labels = {'v_x [km/s]','v_y [km/s]','v_z [km/s]'};
parameter_labels = {'C_R'};

figure();
for k = 1:3
    subplot(3,1,k); hold on
    line = plot(t/3600, +sig3(:,k), lineStyles{1});
    plot(t/3600, -sig3(:,k), lineStyles{1},'Color',line.Color)
    xlabel('Time [hr]')
    ylabel(r_labels{k})
end
sgtitle(filter_type + " Position \pm3\sigma Bounds")
subplot(3,1,1)
lgd = legend({'\pm3\sigma'});
lgd.Units = 'normalized';
lgd.Position = [0.75 0.935 0.25 0.05];

figure();
for k = 1:3
    subplot(3,1,k); hold on
    line = plot(t/3600, +sig3(:,k+3), lineStyles{1});
    plot(t/3600, -sig3(:,k+3), lineStyles{1},'Color',line.Color)
    xlabel('Time [hr]')
    ylabel(v_labels{k})
end
sgtitle(filter_type + " Velocity \pm3\sigma Bounds")
subplot(3,1,1)
lgd = legend({'\pm3\sigma'});
lgd.Units = 'normalized';
lgd.Position = [0.75 0.935 0.25 0.05];

figure();
for k = 1:length(parameter_labels)
    subplot(length(parameter_labels),1,k); hold on
    line = plot(t/3600, +sig3(:,k+6), lineStyles{1});
    plot(t/3600, -sig3(:,k+6), lineStyles{1},'Color',line.Color)
    xlabel('Time [hr]')
    ylabel(parameter_labels{k})
end
sgtitle(filter_type + " Parameter \pm3\sigma Bounds")
subplot(length(parameter_labels),1,1)
lgd = legend({'\pm3\sigma'});
lgd.Units = 'normalized';
lgd.Position = [0.75 0.935 0.25 0.05];

%% Residual plots (pre + post)
if data_type == "range"
    y_labels_time = {'Range Residual [km]'};
    y_labels_QQ = {'Quantiles of Range Residual'};
    sig_res = sqrt(R(1,1)).';
elseif data_type == "range rate"
    y_labels_time = {'Range-Rate Residual [km/s]'};
    y_labels_QQ = {'Quantiles of Range-Rate Residual'};
    sig_res = sqrt(R(2,2)).';
else
    y_labels_time = {'Range Residual [km]','Range-Rate Residual [km/s]'};
    y_labels_QQ = {'Quantiles of Range Residual','Quantiles of Range-Rate Residual'};
    sig_res = sqrt(diag(R)).';
end
sig_res3 = 3.*sig_res;

% Station coloring
stn = station_id(:);
ustn = unique(stn(~isnan(stn)));
ns = numel(ustn);
C = lines(max(ns,1));
function time_scatter_by_station(t_hr, residual_vec, stn_vec, ylabel_str)
    hold on; grid on
    for s = 1:ns
        sid = ustn(s);
        m = (stn_vec == sid) & ~isnan(residual_vec);
        if ~any(m), continue; end
        scatter(t_hr(m), residual_vec(m), 14, C(s,:), 'filled', 'MarkerFaceAlpha', 0.75);
    end
    xlabel('Time [hr]')
    ylabel(ylabel_str)
    if ns > 0
        legend(arrayfun(@(x) "DSS " + string(x), ustn, 'UniformOutput', false), ...
               'Location','northeast')
    end
end

m_meas = size(y,2);   % number of measurement types

% Pre-fit
figure();
for k = 1:m_meas
    % Residual vs time
    nexttile(2*k-1);
    time_scatter_by_station(t/3600, y(:,k), stn, y_labels_time{k});
    hold on

    ax = gca;

    % +/- 3sigma bounds
    ub =  sig_res3(k) * ones(size(t));
    lb = -sig_res3(k) * ones(size(t));
    
    % h3u = plot(t/3600, ub, lineStyles{1}, 'LineWidth',1.5, 'DisplayName','\pm3\sigma');
    % plot(t/3600, lb, lineStyles{1}, 'LineWidth',1.5, 'Color', h3u.Color, ...
    %        'HandleVisibility','off');

    title(filter_type + " Pre-Fit Residuals (" + y_labels_time{k} + ")")

    % Histogram with Gaussian overlay
    nexttile(2*k);
    
    data = y(:,k);
    data = data(~isnan(data));   % remove NaNs if present
    
    histogram(data, 'Normalization','pdf')
    hold on
    
    % Fit normal
    mu = mean(data);
    sigma = std(data);
    
    x = linspace(mu-4*sigma, mu+4*sigma, 100);
    plot(x, normpdf(x,mu,sigma), 'r','LineWidth',1.5)
    
    ylabel('PDF')
    xlabel('Residual')
    title(filter_type + " Pre-Fit Histogram (" + y_labels_time{k} + ")")
    legend('Residuals','Gaussian fit')
    hold off
end

% Post-fit
figure();
for k = 1:m_meas
    % Residual vs time
    nexttile(2*k-1);
    time_scatter_by_station(t/3600, yhat(:,k), stn, y_labels_time{k});
    hold on

    % +/- 3sigma bounds (constant over time)
    ub =  sig_res3(k) * ones(size(t));
    lb = -sig_res3(k) * ones(size(t));

    % h3u = plot(t/3600, ub, lineStyles{1}, 'LineWidth',1.5, 'DisplayName','\pm3\sigma');
    % plot(t/3600, lb, lineStyles{1}, 'LineWidth',1.5, 'Color', h3u.Color, ...
    %            'HandleVisibility','off');

    title(filter_type + " Post-Fit Residuals (" + y_labels_time{k} + ")")

    % Histogram with Gaussian overlay
    nexttile(2*k);
    
    data = yhat(:,k);
    data = data(~isnan(data));   % remove NaNs if present
    
    histogram(data, 'Normalization','pdf')
    hold on
    
    % Fit normal
    mu = mean(data);
    sigma = std(data);
    
    x = linspace(mu-4*sigma, mu+4*sigma, 100);
    plot(x, normpdf(x,mu,sigma), 'r','LineWidth',1.5)
    
    ylabel('PDF')
    xlabel('Residual')
    title(filter_type + " Post-Fit Histogram (" + y_labels_time{k} + ")")
    legend('Residuals','Gaussian fit')
    hold off
end
%% RMS Residuals

rms_prefit = rms(y,'omitnan');
rms_postfit = rms(yhat,'omitnan');

m = size(y,2);

disp("--- " + filter_type + " Pre-fit Residual RMS ---")
for k = 1:m
    disp([y_labels_time{k}, ' RMS = ', num2str(rms_prefit(k))])
end

disp("--- " + filter_type + " Post-fit Residual RMS ---")
for k = 1:m
    disp([y_labels_time{k}, ' RMS = ', num2str(rms_postfit(k))])
end

end
