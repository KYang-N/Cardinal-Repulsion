%% Simulation of the one-layer network with Poisson noise
%% Fisher information

clear
close
SaveFlag = 0;
LoadFlag = 1;
%% Define parameters

if LoadFlag == 0
    seed = 123;
    rng(seed)
    alphaRange = linspace(1e-4,1e-3,4);
    betaRange = linspace(5e-4,5e-3,4);

    MemoryNet.Mode = 'Both';

    MemoryNet.N = 300;

    % Simulation parameters
    DynParams.AddNoise = 1;
    DynParams.dt = 1e-3;
    DynParams.StimTime = 0.5;
    DynParams.RepTime = 1e3;
    DynParams.PFTime = 5;
    DynParams.DecodeTime = [1,2.5,4] + DynParams.StimTime;
    DynParams.Parallel = 2;

%% Dynamics of the full SDE

FIRecord = cell(length(alphaRange),length(betaRange));
p = parpool(DynParams.Parallel);
for ii = 1:length(alphaRange)
    for jj = 1:length(betaRange)
        MemoryNet.alpha = alphaRange(ii);
        MemoryNet.beta = betaRange(jj);
        [FI,~,DynParams] = OneLayerNetworkFI(MemoryNet,DynParams);
        FIRecord{ii,jj} = FI;
    end
end
delete(p);
if SaveFlag
    save('OneLayerFIPatternalphabetaSearch.mat'); %#ok<*UNRCH> 
end
end
%%
if LoadFlag
    clear
    DataDir = '';
    load([DataDir,'OneLayerFIPatternalphabetaSearch.mat']);
end

f = figure;
figure(f)
for ii = 1:4
    for jj = 1:4
        FIMatrix(ii,jj) = min(FIRecord{ii,jj}(end,:)); %#ok<*SAGROW> 
    end
end
imagesc(betaRange,alphaRange,FIMatrix); % Convert to degrees/s
set(gca,'FontSize',10,'Ydir','Normal','TickLength',[0.025,0.01],'TickDir', ...
    'out','LineWidth',0.8,'LooseInset',[0 0 0 0]);
h = colorbar;
h.Ticks = [10 45];
clim([10,45]);
yticks(1e-4:4.5e-4:1e-3);
xticks(5e-4:2.25e-3:5e-3);
xlabel('\beta');
ylabel('\alpha');
axis square
box off;
set(gcf,'Unit','Centimeters','Position',[2,2,7,5]);