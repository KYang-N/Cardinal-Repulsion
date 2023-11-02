function I0 = ExternalInput(NInputSample,Ns,C,epsilon,mu)
NInputSample = NInputSample-1;
dthetas = 2*pi/Ns;
thetas = 0:dthetas:2*pi-dthetas;

dSample = 2*pi/NInputSample;
SampleInput = 0:dSample:2*pi;
NInputSample = NInputSample + 1;

I0 = zeros(NInputSample,Ns);
for kk = 1:NInputSample
    I_ext = C*(1-2*epsilon+2*epsilon*exp(-((thetas-pi)/mu).^2));
    I0(kk,:) = circshift(I_ext,round(((-pi+SampleInput(kk))/dthetas)));
end