function plot_compare_RMS(t,Xtrue,Xhat_list,name_list)
% Inputs:
%   > t: Nx1 time vector [s]
%   > Xtrue : Nx6 true state history [km, km/s]
%   > Xhat_list : cell array of Nx6 estimated state histories (single run each)
%   > name_list : cell array of strings for legend entries
%
% Outputs:
%   > Figures:
%       1) Position component error magnitude vs time: 3 subplots (x,y,z)
%       2) Velocity component error magnitude vs time: 3 subplots (vx,vy,vz)


%% Plot formatting 
set(groot,'defaultFigureColor','w')
set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesLineWidth',1.2)
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesGridAlpha',0.3)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')


%% Compute single-run "RMS" = abs(component error) vs time
M = numel(Xhat_list);
Eabs = cell(1,M); % each Nx6
for i = 1:M
    e = Xhat_list{i} - Xtrue;
    Eabs{i} = abs(e);
end

%% Labels
r_labels = {'x RMS Error [km]','y RMS Error [km]','z RMS Error [km]'};
v_labels = {'v_x RMS Error [km/s]','v_y RMS Error [km/s]','v_z RMS Error [km/s]'};

%% Position RMS plots
figure();
for k = 1:3
    subplot(3,1,k); hold on
    for i = 1:M
        plot(t/3600, Eabs{i}(:,k),'.','MarkerSize',8);
    end
    xlabel('Time [hr]')
    ylabel(r_labels{k})
end
sgtitle('Position RMS Error vs Time')

subplot(3,1,1)
lgd = legend(name_list,'Location','northeast');
lgd.Units = 'normalized';
lgd.Position = [0.72 0.935 0.28 0.05];

%% Velocity RMS plots
figure();
for k = 1:3
    subplot(3,1,k); hold on
    for i = 1:M
        plot(t/3600, Eabs{i}(:,k+3),'.','MarkerSize',8);
    end
    xlabel('Time [hr]')
    ylabel(v_labels{k})
end
sgtitle('Velocity RMS Error vs Time')

subplot(3,1,1)
lgd = legend(name_list,'Location','northeast');
lgd.Units = 'normalized';
lgd.Position = [0.72 0.935 0.28 0.05];

end