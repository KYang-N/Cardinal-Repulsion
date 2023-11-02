% Convert from quantities with respect to length to quantities with respect
% to \theta
function [Dtheta,vtheta] = ConvertTotheta(s,dSample,D,v)
dSample = dSample/2*180/pi; % Convert to degrees
sRight = [s(:,2:end),s(:,2)];
sLeft = [s(:,end-1),s(:,1:end-1)];
ds = sRight-sLeft;
deri = ds/(2*dSample);
vtheta = v./sqrt(sum(deri.^2,1));
Dtheta = D./sum(deri.^2,1);
end