function Amplitude = ConnModulation(Theta,ModuIndex)
% Positive ModuIndex is to strengthen the connectivity at the cardinals
Amplitude = 1+ModuIndex*cos(2*Theta);



