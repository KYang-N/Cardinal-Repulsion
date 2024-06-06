%% Search \alpha and \beta - potential depth and noise coef. index

clear
SaveFlag = 0;
LoadFlag = 1;
%% Define parameters

if LoadFlag == 0
    SensoryNet.N = 300;MemoryNet.N = 300;
    JForward = 0.1;
    nu = 0.17*pi;
    JBackward = 0.25;
    nufdbk = 0.17*pi;

    SensoryNet.Mode = 'Both';
    alphaRange = linspace(0.02,0.08,21);
    betaRange = linspace(0.02,0.12,21);

    SensoryNet.JE = 0.6;
    SensoryNet.JI = 0.35;
    SensoryNet.lambda = 0.36*pi;
    SensoryNet.lambdaI = 1.1*pi;

    % Simulation parameters
    DynParams.dt = 1e-3;
    DynParams.Manifold_tmax = 1.5;
    DynParams.StimTime = 0.5;
    DynParams.NInputSample = 51;
    DynParams.dSample = 2*pi/(DynParams.NInputSample-1);
%% Generate feedforward connectivity

    dthetas = 2*pi/SensoryNet.N;
    thetas = 0:dthetas:2*pi-dthetas;
    SensoryNet.MForward = zeros(MemoryNet.N,SensoryNet.N);
    fForward = 1/(2*pi)*JForward*exp(-((pi-thetas)/nu).^2);
    for i= 1:MemoryNet.N
        SensoryNet.MForward(i,:) = dthetas*circshift(fForward,[0 -round(pi/dthetas)-1+i]);
    end
%% Generate feedback connectivity

    dthetam = 2*pi/MemoryNet.N;
    thetam = 0:dthetam:2*pi-dthetam;
    MemoryNet.MBackward = zeros(SensoryNet.N,MemoryNet.N);
    fBackward = 1/(2*pi)*JBackward*exp(-((pi-thetam)/nufdbk).^2);
    for i= 1:SensoryNet.N
        MemoryNet.MBackward(i,:) = dthetam*circshift(fBackward,[0 -round(pi/dthetam)-1+i]);
    end
%% Evaluate v and D contrast

    tic
    Depth = zeros(length(alphaRange),length(betaRange));
    NCIndex = Depth;
    tic
    for j = 1:length(alphaRange)
        for i = 1:length(betaRange)
            SensoryNet.alpha = alphaRange(j);
            SensoryNet.beta = betaRange(i);
            % Generate recurrent connectivity
            SensoryNet = SensoryNetRecurConn(SensoryNet);
            [vtheta,Dtheta] = Compute_vtheta_Dtheta_new(SensoryNet,MemoryNet,DynParams);
            [Depth(j,i),NCIndex(j,i)] = ComputePotentialDepthNoiseCoefIndex(vtheta,Dtheta,DynParams);
        end
        disp([num2str(j),' alpha value(s) finished.']);
        toc
    end
    toc
     if SaveFlag
        save('TwoLayerDepthNCIndexSearchalphabeta.mat') %#ok<*UNRCH>
    end
end
%%
if LoadFlag
    clear
    DataDir = '';
    load([DataDir,'TwoLayerDepthNCIndexSearchalphabeta.mat'])
end
FigOutDir = '';
f1 = figure;
figure(f1)
imagesc(betaRange,alphaRange,Depth);
h = colorbar;
clim([-40,40]);
h.Ticks = [-40 0 40];
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

f2 = figure;
figure(f2)
imagesc(betaRange,alphaRange,NCIndex);
h = colorbar;
h.TickDirection = 'out';
h.Ticks = [-0.6,0,0.6];
set(h,'LineWidth',0.8);
clim([-0.6 0.6]);
set(gca,'Ydir','normal','FontSize',10,'LineWidth',0.8,'TickDir', ...
    'out','TickLength',[0.02,0.01],'LooseInset',[0 0 0 0]);
ax = gca;
ax.XAxis.Exponent = -2;
ax.YAxis.Exponent = -2;
xlabel('\beta','FontSize',10);
ylabel('\alpha','FontSize',10);
yticks(0.02:0.03:0.08);
xticks(0.02:0.05:0.12);
set(gcf,'Unit','Centimeters','Position',[2,2,5,5]);
box off
axis square