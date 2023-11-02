% Compute ds/d\theta
function deri = ComputeDerivative(s,Parameters)
NInputSample = length(Parameters);
ParametersRight = [Parameters(2:end),Parameters(2)];
ParametersLeft = [Parameters(end-1),Parameters(1:end-1)];
dtheta = ParametersRight - ParametersLeft;
dtheta(dtheta<0) = dtheta(dtheta<0)+2*pi;
dtheta = dtheta/2*180/pi; % Convert to degrees
sRight = [s(:,2:end),s(:,2)];
sLeft = [s(:,end-1),s(:,1:end-1)];
ds = sRight-sLeft;
deri = ds./dtheta;
for i = 1:NInputSample
    deri(:,i) = deri(:,i)/norm(deri(:,i));
end
end