function plot_smoother_diagnostics(t,Xtrue,XhatSmoother,PSmoother)
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
    > State errors vs time with ±3σ covariance bounds for both filter and
      smoother
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
esmth = XhatSmoother - Xtrue;
ersmth = esmth(:,1:3); 
evsmth = esmth(:,4:6); 


% 3-sigma bounds
N = length(t);
sigsmth = zeros(N,6);
for i = 1:N
    sigsmth(i,:) = sqrt(diag(PSmoother(1:6,1:6,i))).';
end
sig3smth = 3*sigsmth;


%% State error plots
r_labels = {'x Error [km]','y Error [km]','z Error [km]'};
v_labels = {'v_x Error [km/s]','v_y Error [km/s]','v_z Error [km/s]'};

figure();
for k = 1:3
    subplot(3,1,k); hold on
    plot(t/3600, ersmth(:,k), markerStyles{1},'MarkerSize',8)
    line = plot(t/3600, +sig3smth(:,k), lineStyles{1});
    plot(t/3600, -sig3smth(:,k), lineStyles{1},'Color',line.Color)
    xlabel('Time [hr]')
    ylabel(r_labels{k})
end
sgtitle("LKF Smoother Position Errors with \pm3\sigma Bounds")
subplot(3,1,1)
lgd = legend({'Error','\pm3\sigma'});
lgd.Units = 'normalized';
lgd.Position = [0.75 0.935 0.25 0.05];

figure();
for k = 1:3
    subplot(3,1,k); hold on
    plot(t/3600, evsmth(:,k), markerStyles{1},'MarkerSize',8)
    line = plot(t/3600, +sig3smth(:,k+3), lineStyles{1});
    plot(t/3600, -sig3smth(:,k+3), lineStyles{1},'Color',line.Color)
    xlabel('Time [hr]')
    ylabel(v_labels{k})
end
sgtitle("LKF Smoother Velocity Errors with \pm3\sigma Bounds")
subplot(3,1,1)
lgd = legend({'Error','\pm3\sigma'});
lgd.Units = 'normalized';
lgd.Position = [0.75 0.935 0.25 0.05];

end
