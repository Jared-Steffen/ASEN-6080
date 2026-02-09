function plot_station_error(t,Xtrue,Xhat,P,filter_type,station_ids)
%{
Inputs:
    >t: Nx1 time vector [s]
    >Xtrue: Nx6 true state history [km, km/s]
    >Xhat: Nx6 estimated state history [km, km/s]
    >P: 6x6xN state covariance history
    >filter_type: string of filter name for plots
    >station_ids: stations ids correspoding to stations in Rs
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

%% State errors
e = Xhat - Xtrue;
ers = e(:,10:end); 

% 3-sigma bounds
N = length(t);
sig = zeros(N,size(ers,2));
for i = 1:N
    sig(i,:) = sqrt(diag(P(10:end,10:end,i))).';
    if ~isreal(sig(i,:))
        disp("Non-real sigma at index i = " + i)
    end
end
sig3 = 3*sig;

%% State error plots
station_labels = {'x_s Error [km]','y_s Error [km]','z_s Error [km]'};

for j = 1:size(ers,2)/3
    figure();
    ks = 3*(j-1) + (1:3);
    for s = 1:3
        k = ks(s);
        subplot(3,1,s); hold on
        plot(t/3600, ers(:,k), markerStyles{1},'MarkerSize',8)
        line = plot(t/3600, +sig3(:,k), lineStyles{1});
        plot(t/3600, -sig3(:,k), lineStyles{1},'Color',line.Color)
        xlabel('Time [hr]')
        ylabel(station_labels{s})
    end
    sgtitle(filter_type + " Station " + num2str(station_ids(j)) + " Position with \pm3\sigma Bounds")
    subplot(3,1,1)
    lgd = legend({'Error','\pm3\sigma'});
    lgd.Units = 'normalized';
    lgd.Position = [0.75 0.935 0.25 0.05];
end

%% RMS 
rms_xsyszs = rms(ers,'omitnan');

for j = 1:size(ers,2)/3
    ks = 3*(j-1) + (1:3);
    disp("--- " + filter_type + " Station " + num2str(station_ids(j)) + " Position RMS Errors [km] ---")
    disp(['RMS x = ', num2str(rms_xsyszs(ks(1)))])
    disp(['RMS y = ', num2str(rms_xsyszs(ks(2)))])
    disp(['RMS z = ', num2str(rms_xsyszs(ks(3)))])
end

end