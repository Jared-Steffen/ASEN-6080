function plot_constant_error(t,Xtrue,Xhat,P,filter_type)
%{
Inputs:
    >t: Nx1 time vector [s]
    >Xtrue: Nx6 true state history [km, km/s]
    >Xhat: Nx6 estimated state history [km, km/s]
    >P: 6x6xN state covariance history
    >filter_type: string of filter name for plots
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
e_const = e(:,7:9); 

% 3-sigma bounds
N = length(t);
sig = zeros(N,3);
for i = 1:N
    sig(i,:) = sqrt(diag(P(7:9,7:9,i))).';
    if ~isreal(sig(i,:))
        disp("Non-real sigma at index i = " + i)
    end
end
sig3 = 3*sig;

%% State error plots
const_labels = {'\mu Error [km]','J2 Error [km]','C_D Error [km]'};

figure();
for k = 1:3
    subplot(3,1,k); hold on
    plot(t/3600, e_const(:,k), markerStyles{1},'MarkerSize',8)
    line = plot(t/3600, +sig3(:,k), lineStyles{1});
    plot(t/3600, -sig3(:,k), lineStyles{1},'Color',line.Color)
    xlabel('Time [hr]')
    ylabel(const_labels{k})
end
sgtitle(filter_type + " Constants with \pm3\sigma Bounds")
subplot(3,1,1)
lgd = legend({'Error','\pm3\sigma'});
lgd.Units = 'normalized';
lgd.Position = [0.75 0.935 0.25 0.05];

%% RMS 
rms_const = rms(e_const,'omitnan');

disp("--- " + filter_type + " Constant RMS Errors [km] ---")
disp(['RMS mu = ', num2str(rms_const(1))])
disp(['RMS J2 = ', num2str(rms_const(2))])
disp(['RMS CD = ', num2str(rms_const(3))])


end