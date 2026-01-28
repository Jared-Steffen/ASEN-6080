function plot_rv_state(t,states,legends,delta_bool)
%{
Inputs:
    >t: time vectors
    >states: cell array of the time histories of the state variables to be
             plotted (corresponding to related time vector)
    >legends: optional list of legend categroies
    >delta_bool: boolean to determine if plots are full and perturbed
Outputs:
    > Plots of the time history evolution of state variables r and v
%}

% Set default plotting settings
set(groot,'defaultFigureColor','w')
set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesLineWidth',1.2)
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesGridAlpha',0.3)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
lineStyles = {'-','--',':','-.'};

if delta_bool
    r_labels = {'\deltax Position [km]','\deltay Position [km]','\deltaz Position [km]'};
    v_labels = {'\deltav_x Velocity [km/s]','\deltav_y Velocity [km/s]','\deltav_z Velocity [km/s]'};
    titles = {'\deltar State History','\deltav State History'};
else
    r_labels = {'x Position [km]','y Position [km]','z Position [km]'};
    v_labels = {'v_x Velocity [km/s]','v_y Velocity [km/s]','v_z Velocity [km/s]'};
    titles = {'r State History','v State History'};
end

figure();
for k = 1:length(states)
    state = states{k};
    ls = lineStyles{k};

    % Plot r components
    r = state(:,1:3);

    subplot(311)
    hold on
    plot(t/3600,r(:,1),'LineStyle',ls)
    xlabel('Time [hr]')
    ylabel(r_labels{1})
    subplot(312)
    hold on
    plot(t/3600,r(:,2),'LineStyle',ls)
    xlabel('Time [hr]')
    ylabel(r_labels{2})
    subplot(313)
    hold on
    plot(t/3600,r(:,3),'LineStyle',ls)
    xlabel('Time [hr]')
    ylabel(r_labels{3})
    sgtitle(titles{1})
end
if ~isempty(legends)
    subplot(311)
    lgd = legend(legends);
    lgd.Units = 'normalized';
    lgd.Position = [0.75 0.935 0.25 0.05];
end

figure();
for k = 1:length(states)
    state = states{k};
    ls = lineStyles{k};

    % Plot v components
    v = state(:,4:6);
    
    subplot(311)
    hold on
    plot(t/3600,v(:,1),'LineStyle',ls)
    xlabel('Time [hr]')
    ylabel(v_labels{1})
    subplot(312)
    hold on
    plot(t/3600,v(:,2),'LineStyle',ls)
    xlabel('Time [hr]')
    ylabel(v_labels{2})
    subplot(313)
    hold on
    plot(t/3600,v(:,3),'LineStyle',ls)
    xlabel('Time [hr]')
    ylabel(v_labels{3})
    sgtitle(titles{2})
end
if ~isempty(legends)
    subplot(311)
    lgd = legend(legends);
    lgd.Units = 'normalized';
    lgd.Position = [0.75 0.935 0.25 0.05];
end

end

