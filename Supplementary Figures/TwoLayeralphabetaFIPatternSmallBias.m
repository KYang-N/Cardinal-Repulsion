%% Simulation of the two-layer network with Poisson noise
%% FI pattern -- different \alpha and \beta

clear
close
DataDir = '';
LoadFlag = 1;
SaveFlag = 0;
%% Define parameters

if LoadFlag == 0
    seed = 123;
    rng(seed)
    load([DataDir,'TwoLayerBiasSDSearchalphabeta.mat']);
    ValidPos = abs(MaxBias)<0.4;
    clearvars -except SaveFlag ValidPos;
    SensoryNet.N = 300;
    alphaRange = linspace(0.02,0.08,20);
    betaRange = linspace(0.02,0.12,20);
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
    DynParams.RepTime = 1e3;
    DynParams.PFTime = 5;
    DynParams.FullDecodeTime = [1,2.5,4] + DynParams.StimTime;
    DynParams.Parallel = 4;

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

    p = parpool(DynParams.Parallel);
    FIRecord = cell(length(alphaRange),length(betaRange));
    for ii = 1:length(alphaRange)
        for jj= 1:length(betaRange)
            if ValidPos(ii,jj) == 1 || ((ii+jj)<=22 && (ii+jj)>=20)
                SensoryNet.alpha = alphaRange(ii);
                SensoryNet.beta = betaRange(jj);
                [FI,~,MemoryNet,DynParams] = TwoLayerNetworkFI(SensoryNet,MemoryNet,DynParams);
                FIRecord{ii,jj} = FI;
            else
                FIRecord{ii,jj} = [];
            end

        end
    end
    delete(p);
    if SaveFlag
        save('TwoLayeralphabetaFIPatternSmallBias.mat'); %#ok<*UNRCH> 
    end
end
%%
if LoadFlag
    clear
    DataDir = '';
    load([DataDir,'TwoLayeralphabetaFIPatternSmallBias.mat']);
end

for ii = 1:20
    for jj = 1:20
        if ~isempty(FIRecord{ii,jj})
            FIMatrix(ii,jj) = min(FIRecord{ii,jj}(end,:)); %#ok<*SAGROW> 
        else
            FIMatrix(ii,jj) = nan;
        end
    end
end

f = figure;
figure(f)
h = heatmap(flipud(FIMatrix),'MissingDataColor',[1 1 1],'ColorLimits',[10 45],'GridVisible','off');
h.CellLabelFormat = '%.1f';
xlabel('\beta');
ylabel('\alpha');
set(gca,'FontSize',10);
set(gcf,'Unit','centimeters','Position',[2,2,18,15]);