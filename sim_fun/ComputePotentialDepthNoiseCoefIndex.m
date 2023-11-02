function [Depth,NCIndex] = ComputePotentialDepthNoiseCoefIndex(vtheta,Dtheta,DynParams)
    InputSample = 0:DynParams.dSample:2*pi;
    MaxInd = round(DynParams.NInputSample/4);
    Potential = EstimatePotential(vtheta,InputSample);
    Depth = Potential(1)-Potential(MaxInd);
    NC = sqrt(2*Dtheta);
    NCIndex = (NC(MaxInd)-NC(1))/(NC(MaxInd)+NC(1));
end