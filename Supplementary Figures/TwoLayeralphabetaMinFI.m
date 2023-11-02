%% Simulation of the two-layer network with Poisson noise
%% FI pattern -- different \alpha and \beta

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
            SensoryNet.alpha = alphaRange(ii);
            SensoryNet.beta = betaRange(jj);

            [FI,~,MemoryNet,DynParams] = TwoLayerNetworkFI(SensoryNet,MemoryNet,DynParams);
            FIRecord{ii,jj} = FI;

        end
    end
    delete(p);
    if SaveFlag
        save('TwoLayeralphabetaFIPattern.mat'); %#ok<*UNRCH> 
    end
end
%%
if LoadFlag
    clear
    DataDir = '';
    load([DataDir,'TwoLayeralphabetaFIPattern.mat']);
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
axis square;
box off;
xlabel('\beta');
ylabel('\alpha');
set(gcf,'Unit','Centimeters','Position',[2,2,7,5]);
yticks(0.02:0.03:0.08);
xticks(0.02:0.05:0.12);
ax = gca;
ax.XAxis.Exponent = -2;
ax.YAxis.Exponent = -2;