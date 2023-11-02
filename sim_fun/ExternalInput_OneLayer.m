function I0 = ExternalInput_OneLayer(NInputSample,Nm)
NInputSample = NInputSample-1;
dthetam = 2*pi/Nm;
thetam = 0:dthetam:2*pi-dthetam;

dSample = 2*pi/NInputSample;
SampleInput = 0:dSample:2*pi;
NInputSample = NInputSample + 1;

I0 = zeros(NInputSample,Nm);
for kk = 1:NInputSample
    Iext = (cos(thetam-pi)+1)/2;
    I0(kk,:) = circshift(Iext,round(((-pi+SampleInput(kk))/dthetam)));
end