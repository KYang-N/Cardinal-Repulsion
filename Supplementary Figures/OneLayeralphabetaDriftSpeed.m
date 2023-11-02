%% Simulation of the one-layer network with Poisson noise - different alpha and beta

clear
close
LoadFlag = 1;
SaveFlag = 0;
%% Define parameters

if LoadFlag == 0
seed = 123;
rng(seed)
alphaRange = linspace(1e-4,1e-3,4);
betaRange = linspace(5e-4,5e-3,4);

MemoryNet.Mode = 'Both';

MemoryNet.N = 300;

% Simulation parameters
DynParams.AddNoise = 0;
DynParams.dt = 1e-3;
DynParams.StimTime = 0.5;
DynParams.PFTime = 5;
DynParams.DecodeTime = [0,4] + DynParams.StimTime;
DynParams.Parallel = 0;
%% Dynamics of the full SDE

DecodedthetaRecord = cell(length(alphaRange),length(betaRange));
for ii = 1:length(alphaRange)
    for jj = 1:length(betaRange)
        MemoryNet.alpha = alphaRange(ii);
        MemoryNet.beta = betaRange(jj);
        [DecodedthetaFull,~,DynParams] = OneLayerDynamics(MemoryNet,DynParams);
        DecodedthetaRecord{ii,jj} = DecodedthetaFull;
    end
end
if SaveFlag 
    save('OneLayeralphabetaDriftSpeed.mat');
end
end
%%
if LoadFlag
    clear
    DataDir = '';
    load([DataDir,'OneLayeralphabetaDriftSpeed.mat']);
end

DriftSpeed = zeros(length(alphaRange),length(betaRange));
for ii = 1:length(alphaRange)
    for jj = 1:length(betaRange)
        thetaDiff = DecodedthetaRecord{ii,jj}(:,:,2)-DecodedthetaRecord{ii,jj}(:,:,1);
        thetaDiff(thetaDiff>pi) = thetaDiff(thetaDiff>pi) - 2*pi;
        thetaDiff(thetaDiff<-pi) = thetaDiff(thetaDiff<-pi) + 2*pi;
        thetaDiff = thetaDiff*180/2/pi; % We map 180 to 2*pi
        DriftSpeed(ii,jj) = max(abs(thetaDiff/4));
    end
end

f1 = figure;
figure(f1)
imagesc(betaRange,alphaRange,DriftSpeed);
h = colorbar;
clim([0,4]);
set(h,'LineWidth',0.8);
h.TickDirection = 'out';
set(gca,'Ydir','normal','FontSize',10,'LineWidth',0.8,'TickDir','out', ...
    'TickLength',[0.02,0.01],'LooseInset',[0 0 0 0]);
xlabel('\beta','FontSize',10);
ylabel('\alpha','FontSize',10);
yticks(1e-4:4.5e-4:1e-3);
xticks(5e-4:2.25e-3:5e-3);
box off
axis square
set(gcf,'Unit','Centimeters','Position',[2,2,5,5]);