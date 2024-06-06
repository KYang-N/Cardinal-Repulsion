%% Simulation of the one-layer network with Poisson noise

clear
close
LoadFlag = 1;
SaveFlag = 0;
%% Define parameters

if LoadFlag == 0
    seed = 123;
    rng(seed)

    MemoryNet.alpha = 4e-4;
    MemoryNet.beta = 2.4e-3;
    MemoryNet.Mode = 'Both'; % Both or EOnly to select if the inhibitory modulation is inlucded

    MemoryNet.N = 300;

    % Simulation parameters
    DynParams.AddNoise = 1;
    DynParams.dt = 1e-3;
    DynParams.StimTime = 0.5;
    DynParams.RepTime = 1e3;
    DynParams.PFTime = 5;
    DynParams.DecodeTime = [1,2.5,4] + DynParams.StimTime;
    DynParams.Parallel = 10;
%% Dynamics of the full SDE

    [DecodedthetaFull,MemoryNet,DynParams] = OneLayerDynamics(MemoryNet,DynParams);

    if SaveFlag
        disp('Saving the data.')
        save(['OneLayerPoissonNoiseResultsalpha',strrep(num2str(MemoryNet.alpha),'.','p'), ...
            'beta',strrep(num2str(MemoryNet.beta),'.','p'),'Mode',num2str(MemoryNet.Mode),'.mat']);
    end
end
%%
if LoadFlag
    clear %#ok<*UNRCH> 
    DataDir = '';
    load([DataDir,'OneLayerPoissonNoiseResultsalpha0p0005beta0p0024ModeBoth.mat']);
end
close
Color = [0.6667 0.6667 0.6667;0.3333 0.3333 0.3333;0 0 0];
SampleInput = 0:DynParams.dSample:2*pi;
f = PlotBiasStdDiscri(DecodedthetaFull,DynParams.DecodeTime,DynParams.StimTime,SampleInput,0,Color);