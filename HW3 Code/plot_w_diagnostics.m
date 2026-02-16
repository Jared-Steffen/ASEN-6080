function plot_w_diagnostics(t, wtrue, Xhat, P, filter_type)
%{
Inputs:
    >t: Nx1 time vector [s]
    >wtrue: Nx3 true J3/DMC acceleration history [km/s^2]
    >Xhat: Nx9 estimated state history (states 7–9 are w_hat)
    >P: 9x9xN covariance history
    >filter_type: string of filter name for plots
Outputs:
    > Acceleration component errors vs time with ±3σ bounds
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

%% Extract estimated w and compute errors
what = Xhat(:,7:9);
ew = what - wtrue;    % Nx3 acceleration error

%% 3-sigma bounds
N = length(t);
sigw = zeros(N,3);

for i = 1:N
    sigw(i,:) = sqrt(diag(P(7:9,7:9,i))).';
    if ~isreal(sigw(i,:))
        disp("Non-real sigma_w at index i = " + i)
    end
end

sigw3 = 3*sigw;

%% Component plots
w_labels = {'w_x Error [km/s^2]', ...
            'w_y Error [km/s^2]', ...
            'w_z Error [km/s^2]'};

figure();
for k = 1:3
    subplot(3,1,k); hold on
    plot(t/3600, ew(:,k), markerStyles{1}, 'MarkerSize', 8)
    line = plot(t/3600, +sigw3(:,k), lineStyles{1});
    plot(t/3600, -sigw3(:,k), lineStyles{1}, 'Color', line.Color)
    xlabel('Time [hr]')
    ylabel(w_labels{k})
end

sgtitle(filter_type + " J3/DMC Acceleration Error with \pm3\sigma Bounds")

subplot(3,1,1)
lgd = legend({'Error','\pm3\sigma'});
lgd.Units = 'normalized';
lgd.Position = [0.75 0.935 0.25 0.05];

end