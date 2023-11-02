%% Simulation of the two-layer network with Poisson noise
%% Bias, SD, and Discrimination Threshold

clear
close
SaveFlag = 0;
LoadFlag = 1;
%% Define parameters

if LoadFlag == 0
    seed = 123;
    rng(seed)

    SensoryNet.N = 300;
    SensoryNet.alpha = 0.03;
    SensoryNet.beta = 0.08;
    SensoryNet.Mode = 'Both';
    MemoryNet.N = 300;

    JForward = .10;
    nu = 0.17*pi;
    JBackward = .25;
    nufdbk = 0.17*pi;

    % Simulation parameters
    DynParams.AddNoise = 1;
    DynParams.dt = 1e-3;
    DynParams.StimTime = 0.5;
    DynParams.RepTime = 3;
    DynParams.PFTime = 5; 
    DynParams.FullDecodeTime = [1,2.5,4] + DynParams.StimTime;
    DynParams.Decoder = 'Aware'; % This sets the preferred orientation to the maxima of the tuning curves at steady states (at PFTime)                       
    DynParams.Parallel = 2;
%% Generate feedforward connectivity

    dthetas = 2*pi/SensoryNet.N;
    thetas = 0:dthetas:2*pi-dthetas;
    SensoryNet.MForward = zeros(MemoryNet.N,SensoryNet.N);
    fForward = 1/(2*pi)*JForward*exp(-((pi-thetas)/nu).^2);
    for i= 1:MemoryNet.N
        SensoryNet.MForward(i,:) = dthetas*circshift(fForward,[0 ...
            -round(pi/dthetas)-1+i]);
    end
%% Generate feedback connectivity

    dthetam = 2*pi/MemoryNet.N;
    thetam = 0:dthetam:2*pi-dthetam;
    MemoryNet.MBackward = zeros(SensoryNet.N,MemoryNet.N);
    fBackward = 1/(2*pi)*JBackward*exp(-((pi-thetam)/nufdbk).^2);
    for i= 1:SensoryNet.N
        MemoryNet.MBackward(i,:) = dthetam*circshift(fBackward,[0 ...
            -round(pi/dthetam)-1+i]);
    end
%% Dynamics of the full SDE

    [DecodedthetaFull,SesoryNet,MemoryNet,DynParams] = FullSDEDynamics(SensoryNet,MemoryNet,DynParams);

    if SaveFlag
        disp('Saving the data.')
        save(['PoissonNoiseResultsalpha',strrep(num2str(SensoryNet.alpha),'.','p'), ...
            'beta',strrep(num2str(SensoryNet.beta),'.','p'),'Mode',SensoryNet.Mode,'.mat']);
    end
end
%%
if LoadFlag
    clear %#ok<*UNRCH> 
    DataDir = '';
    load([DataDir,'PoissonNoiseResultsalpha0p03beta0p08ModeBothUnaware.mat']);
end
close
Color = [0.6667,0.6667,0.6667;0.3333,0.3333,0.3333;0,0,0];
%% Plot bias, SD, and discriminability pattern

SampleInput = 0:DynParams.dSample:2*pi;
f = PlotBiasStdDiscri(DecodedthetaFull,DynParams.FullDecodeTime,DynParams.StimTime,SampleInput,0,Color);