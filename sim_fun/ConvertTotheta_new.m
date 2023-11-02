% Convert from quantities with respect to length to quantities with respect
% to \theta
function [Dtheta,vtheta] = ConvertTotheta_new(s,ThetaLoc,D,v)
ThetaLoc = ThetaLoc/2*180/pi; % Convert to degrees
sRight = [s(:,2:end),s(:,2)]; % note that s(:,1) == s(:,end)
sLeft = [s(:,end-1),s(:,1:end-1)];
ds = sRight-sLeft;
ThetaRight = [ThetaLoc(2:end),ThetaLoc(2)];
ThetaLeft = [ThetaLoc(end-1),ThetaLoc(1:end-1)];
dtheta = ThetaRight-ThetaLeft;
dtheta(dtheta<0) = dtheta(dtheta<0) + 180;
deri = ds./dtheta;
vtheta = v./sqrt(sum(deri.^2,1));
Dtheta = D./sum(deri.^2,1);
end