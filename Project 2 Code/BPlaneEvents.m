function [value,isterminal,direction] = BPlaneEvents(t,s,RSOI_gate)

r = s(1:3);
value = norm(r) - RSOI_gate; % event when value = 0
isterminal = 1; % stop integration at the event
direction = 0;

end

