function plot_filter_diagnostics(t,Xtrue,Xhat,P,yhat,filter_type)
%{
Inputs:
    >t: Nx1 time vector [s]
    >Xtrue: Nx6 true state history [km, km/s]
    >Xhat: Nx6 estimated state history [km, km/s]
    >P: 6x6xN state covariance history
    >yhat: cell(N,1) OR numeric (Nxm) post-fit residuals
    >R: measurement noise covariance matrix
    >filter_type: string of filter name for plots
Outputs:
    > State errors vs time with ±3σ covariance bounds
    > Post-fit measurement residuals with noise comparison
%}

% Plot formatting
set(groot,'defaultFigureColor','w')
set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesLineWidth',1.2)
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesGridAlpha',0.3)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
lineStyles = {'-','--',':','-.'};
markerStyles = {'.','o','x'};

% State errors
e  = Xhat - Xtrue;     % Nx6
er = e(:,1:3);        % position error
ev = e(:,4:6);        % velocity error

er3 = sqrt(sum(er.^2,2));  % 3D position error
ev3 = sqrt(sum(ev.^2,2));  % 3D velocity error

% 3 sigma bound
N = length(t);
sig = zeros(N,6);
for i = 1:N
    sig(i,:) = sqrt(diag(P(:,:,i))).';
end
sig3 = 3*sig;

% State error plots
r_labels = {'x Error [km]','y Error [km]','z Error [km]'};
v_labels = {'v_x Error [km/s]','v_y Error [km/s]','v_z Error [km/s]'};

figure();
for k = 1:3
    subplot(3,1,k); hold on
    plot(t/3600, er(:,k), markerStyles{1})
    line = plot(t/3600, +sig3(:,k), lineStyles{2});
    plot(t/3600, -sig3(:,k), lineStyles{2},'Color',line.Color)
    xlabel('Time [hr]')
    ylabel(r_labels{k})
end
sgtitle(filter_type + " Position Error with \pm3\sigma Bounds")
subplot(3,1,1)
lgd = legend({'Error','\pm3\sigma'});
lgd.Units = 'normalized';
lgd.Position = [0.75 0.935 0.25 0.05];

figure();
for k = 1:3
    subplot(3,1,k); hold on
    plot(t/3600, ev(:,k), markerStyles{1});
    line = plot(t/3600, +sig3(:,k+3), lineStyles{2});
    plot(t/3600, -sig3(:,k+3), lineStyles{2},'Color',line.Color)
    xlabel('Time [hr]')
    ylabel(v_labels{k})
end
sgtitle(filter_type + " Velocity Error with \pm3\sigma Bounds")
subplot(3,1,1)
lgd = legend({'Error','\pm3\sigma'});
lgd.Units = 'normalized';
lgd.Position = [0.75 0.935 0.25 0.05];

% Residual plot labels
y_labels = {'Range [km]','Range Rate [km/s]'};
y_labels_QQ = {'Quantiles of Range','Quantiles of Range Rate'};
titles_QQ = {filter_type + " QQ Plot of Range Residuals vs Standard Normal",filter_type + " QQ Plot of Range Rate Residuals vs Standard Normal"};

% Post fit residuals
figure()
for k = 1:2
    subplot(2,1,k); hold on
    plot(t/3600, yhat(:,k), markerStyles{1})
    xlabel('Time [hr]')
    ylabel(y_labels{k})
end
sgtitle(filter_type + " Measurement Post-Fit Residuals")

% Post fit residuals QQ Plot
figure()
for k = 1:2
    subplot(2,1,k); hold on
    qqplot(yhat(:,k))
    ylabel(y_labels_QQ{k})
    title(titles_QQ{k})
end

% RMS
rms_xyz = rms(er);
rms_vxvyvz = rms(ev);
rms_r = rms(er3);
rms_v = rms(ev3);
rms_res = rms(yhat,'omitnan');

disp("--- " + filter_type + " Position RMS Errors [km] ---")
disp(['RMS x = ', num2str(rms_xyz(1))])
disp(['RMS y = ', num2str(rms_xyz(2))])
disp(['RMS z = ', num2str(rms_xyz(3))])
disp(['3D Position RMS = ', num2str(rms_r)])
disp("--- " + filter_type + " Velocity RMS Errors [km/s] ---")
disp(['RMS vx = ', num2str(rms_vxvyvz(1))])
disp(['RMS vy = ', num2str(rms_vxvyvz(2))])
disp(['RMS vz = ', num2str(rms_vxvyvz(3))])
disp(['3D Velocity RMS = ', num2str(rms_v)])
disp("--- " + filter_type + " Post-fit Residual RMS [km]/[km/s] ---")
disp(['Range Residual RMS = ', num2str(rms_res(1))])
disp(['Range Rate Residual RMS = ', num2str(rms_res(2))])


end

