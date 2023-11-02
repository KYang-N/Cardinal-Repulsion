%% Simulation of the two-layer network with Poisson noise
%% Bias, SD -- different \alpha and \beta

clear
close
LoadFlag = 1;
SaveFlag = 0;
%% Define parameters

if LoadFlag == 0
    seed = 123;
    rng(seed)

    SensoryNet.N = 300;
    alphaRange = linspace(0.02,0.08,4);
    betaRange = linspace(0.02,0.12,4);
    SensoryNet.Mode = 'Both';
    MemoryNet.N = 300;

    JForward = .10;
    nu = 0.17*pi;
    JBackward = .25;
    nufdbk = 0.17*pi;

    % Simulation parameters
    DynParams.AddNoise = 0;
    DynParams.dt = 1e-3;
    DynParams.StimTime = 0.5;
    DynParams.PFTime = 5;
    DynParams.FullDecodeTime = [0 4] + DynParams.StimTime;
    DynParams.Decoder = 'Aware';
    DynParams.Parallel = 0;
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

    DecodedthetaRecord = cell(length(alphaRange),length(betaRange));
    for ii = 1:length(alphaRange)
        for jj= 1:length(betaRange)
            SensoryNet.alpha = alphaRange(ii);
            SensoryNet.beta = betaRange(jj);
            [DecodedthetaFull,~,~,DynParams] = FullSDEDynamics(SensoryNet,MemoryNet,DynParams);
            DecodedthetaRecord{ii,jj} = DecodedthetaFull;
        end
    end
    if SaveFlag
        save('TwoLayeralphabetaDriftSpeed.mat'); %#ok<*UNRCH> 
    end
end
%%
if LoadFlag
    clear
    DataDir = '';
    load([DataDir,'TwoLayeralphabetaDriftSpeed.mat']);
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
clim([0 4]);
set(h,'LineWidth',0.8);
h.TickDirection = 'out';
set(gca,'Ydir','normal','FontSize',10,'LineWidth',0.8,'TickDir','out', ...
    'TickLength',[0.02,0.01],'LooseInset',[0 0 0 0]);
xlabel('\beta','FontSize',10);
ylabel('\alpha','FontSize',10);
yticks(0.02:0.03:0.08);
xticks(0.02:0.05:0.12);
ax = gca;
ax.XAxis.Exponent = -2;
ax.YAxis.Exponent = -2;
box off
axis square
set(gcf,'Unit','Centimeters','Position',[2,2,5,5]);