function plot_covariance_trace(t,P,filter_type)
% plot_covariance_trace
%{
Inputs:
    >t: Nx1 time vector [s]
    >P: 6x6xN state covariance history
    >filter_type: string label (e.g. "Batch", "CKF", "EKF")
Outputs:
    > Trace of position and velocity covariance vs time
%}

%% Preallocate
N = length(t);
trace_pos = zeros(N,1);
trace_vel = zeros(N,1);

%% Compute traces
for k = 1:N
    Pk = P(:,:,k);
    trace_pos(k) = trace(Pk(1:3,1:3));   % position covariance trace
    trace_vel(k) = trace(Pk(4:6,4:6));   % velocity covariance trace
end

%% Plot
figure();

subplot(2,1,1); hold on; grid on
plot(t/3600, trace_pos,'-o', 'LineWidth', 2)
xlabel('Time [hr]')
ylabel('tr(P_{rr}) [km^2]')
title(filter_type + " Position Covariance Trace")

subplot(2,1,2); hold on; grid on
plot(t/3600, trace_vel,'-o', 'LineWidth', 2)
xlabel('Time [hr]')
ylabel('tr(P_{vv}) [(km/s)^2]')
title(filter_type + " Velocity Covariance Trace")

end
