function [t,flattenedObsData] = flattenObsData(obsData,station_ids)
% This functions flatterns observation data that has no overlaps into a
% matrix of [time since epoch, station_id, range, range-rate]


% Get time
t = obsData(:,1);
data = obsData(:,2:end);

% Find NaNs
mask_NaN = ~isnan(data);

% Build new matrix
for i = 1:size(data,1)
    flattenedObsData(i,1) = t(i);
    if mask_NaN(i,1) == true && mask_NaN(i,4) == true
        flattenedObsData(i,2) = station_ids(1);
        flattenedObsData(i,3:4) = data(i,mask_NaN(i,:));
    elseif mask_NaN(i,2) == true && mask_NaN(i,5) == true
        flattenedObsData(i,2) = station_ids(2);
        flattenedObsData(i,3:4) = data(i,mask_NaN(i,:));
    elseif mask_NaN(i,3) == true && mask_NaN(i,6) == true
        flattenedObsData(i,2) = station_ids(3);
        flattenedObsData(i,3:4) = data(i,mask_NaN(i,:));
    else
        error('Error: column mismatch')
    end

end

end

