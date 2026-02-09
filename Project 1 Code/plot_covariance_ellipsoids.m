function plot_covariance_ellipsoids(P, filter_type)
%{
Inputs:
    >P: 6x6xN covariance history OR cell array of covariance histories
    >filter_type: string label OR cell/string array of labels
Outputs:
    > 3sigma 2D position covariance ellipses (x-y, x-z, y-z) at final time
    > Allows multiple filters to be overlaid for comparison
%}

% Allow single or multiple inputs
if ~iscell(P)
    P = {P};
end

if ~iscell(filter_type)
    filter_type = cellstr(filter_type);
end

% If one label but multiple covariances, auto-number
if numel(filter_type) == 1 && numel(P) > 1
    base = string(filter_type{1});
    filter_type = cellstr(base + " " + string(1:numel(P)));
end

nsig = 3;
npts = 200;
th = linspace(0,2*pi,npts);

% Extract final-time position covariances
M = numel(P);
Pr = zeros(3,3,M);

for i = 1:M
    PF = P{i}(:,:,end);
    Pr(:,:,i) = PF(1:3,1:3);  % position block
end

% Coordinate pairs: x-y, x-z, y-z
pairs = {[1 2],'x-y','x [km]','y [km]'; ...
         [1 3],'x-z','x [km]','z [km]'; ...
         [2 3],'y-z','y [km]','z [km]'};

figure();
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for k = 1:3
    ij   = pairs{k,1};
    name = pairs{k,2};
    xl   = pairs{k,3};
    yl   = pairs{k,4};

    nexttile; hold on; grid on; axis equal;
    title(name)
    xlabel(xl); ylabel(yl);

    for i = 1:M
        % 2x2 covariance (no clamping)
        C2 = Pr(ij,ij,i);

        % Eigen-decomposition ellipse
        [V,D] = eig(C2);
        A = V * diag(sqrt(diag(D))) * nsig;

        xy = A * [cos(th); sin(th)];
        plot(xy(1,:), xy(2,:), 'LineWidth', 1.5);
    end

    legend(string(filter_type), 'Location','northeast');
    sgtitle('Covariance Ellipsoids')

end

end
