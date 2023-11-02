function Discrimination = ComputeDT(Bias,SqrtVar,SampleInput)
% Bias, Std, and SampleInput are row vectors
% Bias in degree, Std in degree, SampleInput in rad
Dalpha = 1;
SampleInput = SampleInput/pi*180;
dSample = SampleInput(2) - SampleInput(1);
if SampleInput(1) == SampleInput(end)-360
    bprime = (Bias(2:end)-Bias(1:end-1))/dSample;
    Discrimination = Dalpha*SqrtVar(1:end-1)./(1+bprime);
    Discrimination(end+1) = Discrimination(1);
else
    error('The beginning and ending orientations must have a difference of 2*pi.');
end