%% Simulation of the one-layer network with Poisson noise
%% Fisher information -- small bias

clear
DataDir = '';
SaveFlag = 0;
LoadFlag = 1;
%% Define parameters

if LoadFlag == 0
    load([DataDir,'OneLayerBiasSDSearchalphabeta.mat']);
    ValidPos = abs(MaxBias)<0.2;
    clearvars -except SaveFlag ValidPos;
    seed = 123;
    rng(seed)
    alphaRange = linspace(1e-4,1e-3,20);
    betaRange = linspace(5e-4,5e-3,20);
    load()
    MemoryNet.Mode = 'Both';

    MemoryNet.N = 300;

    % Simulation parameters
    DynParams.AddNoise = 1;
    DynParams.dt = 1e-3;
    DynParams.StimTime = 0.5;
    DynParams.RepTime = 1e3;
    DynParams.PFTime = 5;
    DynParams.DecodeTime = [1,2.5,4] + DynParams.StimTime;
    DynParams.Parallel = 6;

%% Dynamics of the full SDE

    FIRecord = cell(length(alphaRange),length(betaRange));
    p = parpool(DynParams.Parallel);
    for ii = 1:length(alphaRange)
        for jj = 1:length(betaRange)
            if ValidPos(ii,jj) == 1 || ((ii+jj)<=22 && (ii+jj)>=20)
                MemoryNet.alpha = alphaRange(ii);
                MemoryNet.beta = betaRange(jj);
                [FI,~,DynParams] = OneLayerNetworkFI(MemoryNet,DynParams);
                FIRecord{ii,jj} = FI;
            else
                FIRecord{ii,jj} = [];
            end
        end
    end
    delete(p);
    if SaveFlag
        save('OneLayerFIPatternalphabetaSearchSmallBias.mat'); %#ok<*UNRCH> 
    end
end
%%
if LoadFlag
    clear
    DataDir = '';
    load([DataDir,'OneLayerFIPatternalphabetaSearchSmallBias.mat']);
end

f2 = figure;
figure(f2)
for ii = 1:20
    for jj = 1:20
        if ~isempty(FIRecord{ii,jj})
            FIMatrix(ii,jj) = min(FIRecord{ii,jj}(end,:)); %#ok<*SAGROW> 
        else
            FIMatrix(ii,jj) = nan;
        end
    end
end

f3 = figure;
figure(f3)
h = heatmap(flipud(FIMatrix),'MissingDataColor',[1 1 1],'ColorLimits',[10 45],'GridVisible','off');
h.CellLabelFormat = '%.1f';
xlabel('\beta');
ylabel('\alpha');
set(gca,'FontSize',10);
set(gcf,'Unit','centimeters','Position',[2,2,18,15]);