function plot_consider_analysis(t,Xtrue,Xhat,P,Pxx)
%{
Inputs:
    >t: Nx1 time vector [s]
    >Xtrue: Nx6 true state history [km, km/s]
    >Xhat: Nx6 estimated state history [km, km/s]
    >P: 6x6xN state covariance history
    >P: 6x6xN state consider covariance history

Outputs:
    > State errors vs time with ±2σ covariance bounds
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
er = e(:,1:3); 
ev = e(:,4:6); 
er3 = sqrt(sum(er.^2,2));
ev3 = sqrt(sum(ev.^2,2));

% 2-sigma bounds w/o consider
N = length(t);
if ~isempty(P)
    sig = zeros(N,6);
    for i = 1:N
        sig(i,:) = sqrt(diag(P(1:6,1:6,i))).';
        if ~isreal(sig(i,:))
            disp("Non-real sigma at index i = " + i)
        end
    end
else
    sig = NaN(N,6);
end
sig2 = 2*sig;

% 2-sigma bounds w/ consider
if ~isempty(Pxx)
    sigc = zeros(N,6);
    for i = 1:N
        sigc(i,:) = sqrt(diag(Pxx(1:6,1:6,i))).';
        if ~isreal(sigc(i,:))
            disp("Non-real sigma for consider at index i = " + i)
        end
    end
else
    sigc = NaN(N,6);
end
sig2c = 2*sigc;



%% State error plots
r_labels = {'x Error [km]','y Error [km]','z Error [km]'};
v_labels = {'v_x Error [km/s]','v_y Error [km/s]','v_z Error [km/s]'};

figure();
for k = 1:3
    subplot(3,1,k); hold on
    plot(t/3600, er(:,k), markerStyles{1},'MarkerSize',8)
    line = plot(t/3600, +sig2(:,k), lineStyles{1});
    plot(t/3600, -sig2(:,k), lineStyles{1},'Color',line.Color)
    line = plot(t/3600, +sig2c(:,k), lineStyles{1});
    plot(t/3600, -sig2c(:,k), lineStyles{1},'Color',line.Color)
    xlabel('Time [hr]')
    ylabel(r_labels{k})
end
sgtitle("Position Error with \pm2\sigma Bounds")
subplot(3,1,1)
lgd = legend('Error','\pm2\sigma: no consider covariance analaysis','','\pm2\sigma: consider covariance analaysis','');
lgd.Units = 'normalized';
lgd.Position = [0.75 0.935 0.25 0.05];

figure();
for k = 1:3
    subplot(3,1,k); hold on
    plot(t/3600, ev(:,k), markerStyles{1},'MarkerSize',8);
    line = plot(t/3600, +sig2(:,k+3), lineStyles{1});
    plot(t/3600, -sig2(:,k+3), lineStyles{1},'Color',line.Color)
    line = plot(t/3600, +sig2c(:,k+3), lineStyles{1});
    plot(t/3600, -sig2c(:,k+3), lineStyles{1},'Color',line.Color)
    xlabel('Time [hr]')
    ylabel(v_labels{k})
end
sgtitle("Velocity Error with \pm2\sigma Bounds")
subplot(3,1,1)
lgd = legend('Error','\pm2\sigma: no consider covariance analaysis','','\pm2\sigma: consider covariance analaysis','');
lgd.Units = 'normalized';
lgd.Position = [0.75 0.935 0.25 0.05];


%% RMS 
rms_xyz = rms(er,'omitnan');
rms_vxvyvz = rms(ev,'omitnan');
rms_r = rms(er3,'omitnan');
rms_v = rms(ev3,'omitnan');

disp("--- Position RMS Errors [km] ---")
disp(['RMS x = ', num2str(rms_xyz(1))])
disp(['RMS y = ', num2str(rms_xyz(2))])
disp(['RMS z = ', num2str(rms_xyz(3))])
disp(['3D Position RMS = ', num2str(rms_r)])

disp("--- Velocity RMS Errors [km/s] ---")
disp(['RMS vx = ', num2str(rms_vxvyvz(1))])
disp(['RMS vy = ', num2str(rms_vxvyvz(2))])
disp(['RMS vz = ', num2str(rms_vxvyvz(3))])
disp(['3D Velocity RMS = ', num2str(rms_v)])

%% Percent outside 2-sigma bounds

state_names = {'x','y','z','vx','vy','vz'};

fprintf('\n--- Per-state percent outside 2-sigma ---\n');
for k = 1:6
    pct_no_consider = 100 * sum(abs(e(:,k)) > sig2(:,k)) / N;
    pct_consider    = 100 * sum(abs(e(:,k)) > sig2c(:,k)) / N;

    fprintf('%3s: no consider = %6.2f%% | consider = %6.2f%%\n', ...
        state_names{k}, pct_no_consider, pct_consider);
end

end
