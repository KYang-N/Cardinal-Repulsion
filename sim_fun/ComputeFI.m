function FI = ComputeFI(SampleInput,TCMean,TCVar)
% TCMean is a N x NInputSample matrix of tuning curves. The returned value
% is a 1 x NInputSample array of total Fisher information.
dSample = SampleInput(2)-SampleInput(1); 
dSample = dSample/2/pi*180; % Convert to degrees
if SampleInput(1) == SampleInput(end)-2*pi
    TCVarMean = (TCVar(:,1:end-1) + TCVar(:,2:end))/2;
    Fisher = ((TCMean(:,2:end)-TCMean(:,1:end-1))/dSample).^2./TCVarMean;
    FI = sum(Fisher,1);
    FI(end+1) = FI(1);
else
    error('The beginning and ending orientations must have a difference of 2*pi.')
end