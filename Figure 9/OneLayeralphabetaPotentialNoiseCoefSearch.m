%% Search \alpha and \beta - potential depth and noise coef. index for the one-layer model

clear
SaveFlag = 0;
LoadFlag = 1;
%% Define parameters

if LoadFlag == 0
    MemoryNet.tau = 0.01;
    MemoryNet.N = 300;
    MemoryNet.Mode = 'Both';
    MemoryNet.NEM = 1.5; MemoryNet.thM = 0.1; MemoryNet.sigM = 6.6; MemoryNet.maxf = 100;
    MemoryNet.q = @(x) MemoryNet.maxf*(x-MemoryNet.thM).^MemoryNet.NEM./( ...
        MemoryNet.sigM^MemoryNet.NEM+(x-MemoryNet.thM).^MemoryNet.NEM).*(x>MemoryNet.thM);

    MemoryNet.JE = 1; MemoryNet.JI = 0.17;
    MemoryNet.lambdaM = 0.2*pi;
    MemoryNet.IEc = 0.6*ones(MemoryNet.N,1);

    alphaRange = linspace(1e-4,1e-3,21);
    betaRange = linspace(5e-4,5e-3,21);

    % Simulation parameters
    DynParams.dt = 1e-3;
    DynParams.Manifold_tmax = 1.5;
    DynParams.StimTime = .5;
    DynParams.NInputSample = round(MemoryNet.N/6)+1;
    DynParams.dSample = 2*pi/(DynParams.NInputSample-1);
%% Evaluate vmax and D contrast

    tic
    Depth = zeros(length(alphaRange),length(betaRange));
    NCIndex = Depth;
    for j = 1:length(alphaRange)
        for i = 1:length(betaRange)
            MemoryNet.alpha = alphaRange(j);
            MemoryNet.beta = betaRange(i);
            % Generate recurrent connectivity
            MemoryNet = OneLayerRecurConn(MemoryNet);
            [vtheta,Dtheta] = Compute_vtheta_Dtheta_OneLayer(MemoryNet,DynParams);
            [Depth(j,i),NCIndex(j,i)] = ComputePotentialDepthNoiseCoefIndex(vtheta,Dtheta,DynParams);
        end
        disp([num2str(j),'alpha value(s) finished.']);
        toc
    end
    
     if SaveFlag
        save('OneLayerDepthNCIndexSearchalphabeta.mat') %#ok<*UNRCH>
    end
end
%%
if LoadFlag
    clear
    DataDir = '';
    load([DataDir,'OneLayerDepthNCIndexSearchalphabeta.mat']);
end

f1 = figure;
figure(f1)
imagesc(betaRange,alphaRange,Depth); 

h = colorbar;
h.Ticks = [-60 0 120];
clim([-60,120]);
set(h,'LineWidth',0.8);
h.TickDirection = 'out';
set(gca,'Ydir','normal','FontSize',10,'LineWidth',0.8,'TickDir','out', ...
    'TickLength',[0.02,0.01],'LooseInset',[0 0 0 0]);
xlabel('\beta','FontSize',10);
ylabel('\alpha','FontSize',10);
yticks(1e-4:4.5e-4:1e-3);
xticks(5e-4:2.25e-3:5e-3);
ax = gca;
ax.XAxis.Exponent = -3;
ax.YAxis.Exponent = -4;
box off
axis square
set(gcf,'Unit','Centimeters','Position',[2,2,5,5]);


f2 = figure;
figure(f2)
imagesc(betaRange,alphaRange,NCIndex);
h = colorbar;
clim([-0.1 0.3])
h.TickDirection = 'out';
h.Ticks = [-0.1,0,0.3];
set(h,'LineWidth',0.8);
% clim([-0.1 0.4]);
set(gca,'Ydir','normal','FontSize',10,'LineWidth',0.8,'TickDir', ...
    'out','TickLength',[0.02,0.01],'LooseInset',[0 0 0 0]);
ax = gca;
ax.XAxis.Exponent = -3;
ax.YAxis.Exponent = -4;
xlabel('\beta','FontSize',10);
ylabel('\alpha','FontSize',10);
yticks(1e-4:4.5e-4:1e-3);
xticks(5e-4:2.25e-3:5e-3);
set(gcf,'Unit','Centimeters','Position',[2,2,5,5]);
box off
axis square
