function plot_measurements(t,measurements,station_ids)
%{
Inputs:
    >t: time vector
    >measurements: cell array of the time histories of the measurements to
                   be plotted [station id, rho, rho_dot, elevation]
    >station_ids: station ids for legend
Outputs:
    > Plots of the time history evolution of rho, rho_dot, and elevation
      measurements corresponding to specific stations
%}

% Set default plotting settings
set(groot,'defaultFigureColor','w')
set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesLineWidth',1.2)
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesGridAlpha',0.3)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

% Set of colors to assign station ids to (up to 12)
colors = [ ...
    0.0000, 0.4470, 0.7410;  % blue
    0.8500, 0.3250, 0.0980;  % orange
    0.9290, 0.6940, 0.1250;  % yellow
    0.4940, 0.1840, 0.5560;  % purple
    0.4660, 0.6740, 0.1880;  % green
    0.3010, 0.7450, 0.9330;  % cyan
    0.6350, 0.0780, 0.1840;  % dark red
    0.9059, 0.1608, 0.5412;  % magenta
    0.1059, 0.6196, 0.4667;  % teal
    0.9020, 0.6706, 0.0078;  % gold
    0.6500, 0.6500, 0.6500;  % gray
    0.4000, 0.4000, 0.8000;  % slate blue
];

% Generate legend
station_legend = cell(length(station_ids),1);
for i = 1:length(station_ids)
    station_legend{i} = "Station ID: "+ station_ids(i);
end

Ns = size(colors,1);
T  = numel(t);

% Initialize
station_id_all = [];
time_all = [];
rho_all = [];
rhodot_all = [];
el_all = [];

for i = 1:T
    M = measurements{i};
    if isempty(M) 
        continue; 
    end
    n = size(M,1);
    station_id_all = [station_id_all; M(:,1)];
    time_all = [time_all; repmat(t(i), n, 1)];
    rho_all = [rho_all; M(:,2)];
    rhodot_all = [rhodot_all; M(:,3)];
    el_all = [el_all; M(:,4)];
end

% Range and range rate figure
figure();
subplot(211) 
hold on
for s = 1:Ns
    idx = (station_id_all == s);
    if any(idx)
        plot(time_all(idx)/3600, rho_all(idx), 'o', ...
             'Color', colors(s,:), 'MarkerSize', 4);
    end
end
ylabel('Range [km]')
xlabel('Time [hr]')
title('Range')

subplot(212)
hold on
for s = 1:Ns
    idx = (station_id_all == s);
    if any(idx)
        plot(time_all(idx)/3600, rhodot_all(idx), 'o', ...
             'Color', colors(s,:), 'MarkerSize', 4);
    end
end
ylabel('Range Rate [km/s]')
xlabel('Time [hr]')
title('Range Rate')
lgd = legend(station_legend);
lgd.Units = 'normalized';
lgd.Position = [0.75 0.935 0.25 0.05];

% Elevation figure
figure();
hold on
for s = 1:Ns
    idx = (station_id_all == s);
    if any(idx)
        plot(time_all(idx)/3600, el_all(idx), 'o', ...
             'Color', colors(s,:), 'MarkerSize', 4);
    end
end
ylabel('Elevation [degrees]')
xlabel('Time [hr]')
title('Elevation')
legend(station_legend)

end

