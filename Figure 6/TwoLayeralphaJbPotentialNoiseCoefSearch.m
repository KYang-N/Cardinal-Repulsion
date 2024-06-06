%% Search \alpha and Jb - drift velocity and coefficient of diffusion of \theta

clear
SaveFlag = 0;
LoadFlag = 1;
%% Define parameters

if LoadFlag == 0
    SensoryNet.tau = 0.01;
    SensoryNet.N = 300;

    JForward = 0.1;
    nu = 0.17*pi;
    nufdbk = 0.17*pi;


    MemoryNet.N = 300;
    alphaRange = linspace(0,0.06,21);

    SensoryNet.Mode = 'EOnly';
    SensoryNet.JI = 0.35;
    SensoryNet.JE = 0.6;
    SensoryNet.lambda = 0.36*pi;
    SensoryNet.lambdaI = 1.1*pi;
    JbRange = linspace(0.15,0.3,21);
    
    % External input
    DynParams.C = 4;
    DynParams.epsilon = 0.2;
    DynParams.mu = 0.3*pi;

% Simulation parameters
DynParams.NInputSample = round(SensoryNet.N/6)+1;
DynParams.dSample = 2*pi/(DynParams.NInputSample-1);
DynParams.AddNoise = 0;
DynParams.dt = 1e-3;
DynParams.StimTime = 0.5;
DynParams.Manifold_tmax = 1.5;
DynParams.PFTime = 1.5;
%% Generate feedforward connectivity

    dthetas = 2*pi/SensoryNet.N;
    thetas = 0:dthetas:2*pi-dthetas;
    SensoryNet.MForward = zeros(MemoryNet.N,SensoryNet.N);
    fForward = 1/(2*pi)*JForward*exp(-((pi-thetas)/nu).^2);
    for i= 1:MemoryNet.N
        SensoryNet.MForward(i,:) = dthetas*circshift(fForward,[0 ...
            -round(pi/dthetas)-1+i]);
    end
%% Evaluate vmax and D contrast

    tic
    PotentialDepth = zeros(length(alphaRange),length(JbRange));
    NCIndex = PotentialDepth;
    for j = 1:length(alphaRange)
        for i = 1:length(JbRange)
            SensoryNet.alpha = alphaRange(j);
            JBackward = JbRange(i);
            % Generate feedback connectivity
            dthetam = 2*pi/MemoryNet.N;
            thetam = 0:dthetam:2*pi-dthetam;
            MemoryNet.MBackward = zeros(SensoryNet.N,MemoryNet.N);
            fBackward = 1/(2*pi)*JBackward*exp(-((pi-thetam)/nufdbk).^2);
            for ii = 1:SensoryNet.N
                MemoryNet.MBackward(ii,:) = dthetam*circshift(fBackward,[0 -round(pi/dthetam)-1+ii]);
            end
    
            % Generate recurrent connectivity
            SensoryNet = SensoryNetRecurConn(SensoryNet);
            [vtheta,Dtheta] = Compute_vtheta_Dtheta_new(SensoryNet,MemoryNet,DynParams);
            [PotentialDepth(j,i),NCIndex(j,i)] = ComputePotentialDepthNoiseCoefIndex( ...
                vtheta,Dtheta,DynParams);
        end
    end
    toc
     if SaveFlag
        save('PotentialNCIndexSearchalphaJb.mat') %#ok<*UNRCH>
    end
end
%%
if LoadFlag
    clear 
    DataDir = '';
    load([DataDir,'PotentialNCIndexSearchalphaJb.mat'])
end
%%
f1 = figure;
figure(f1)
subplot(1,2,1)
imagesc(JbRange,alphaRange,PotentialDepth); % Convert to degrees/s
h = colorbar;
h.TickDirection = 'out';
set(h,'LineWidth',0.8);
h.Ticks = [0 30 60];
xticks([0.15 0.2 0.3]);
yticks([0 0.03 0.06])
clim([0 60])
set(gca,'Ydir','normal','FontSize',10,'LineWidth',0.8,'TickDir','out','TickLength',[0.025,0.01]);
xlabel('$J_b$','FontSize',10,'Interpreter','latex');
ylabel('$\alpha$','FontSize',10,'Interpreter','latex');
box off
axis square
box off

subplot(1,2,2)
imagesc(JbRange,alphaRange,NCIndex);
h = colorbar;
h.TickDirection = 'out';
h.Ticks = [0 0.25 0.5];
set(h,'LineWidth',0.8);
clim([0 0.5])
xticks([0.15 0.2 0.3]);
yticks([0 0.03 0.06])
set(gca,'Ydir','normal','FontSize',10,'LineWidth',0.8,'TickDir','out','TickLength',[0.025,0.01]);
xlabel('$J_b$','FontSize',10,'Interpreter','latex');
% ylabel('$\alpha$','FontSize',10,'Interpreter','latex');

box off
axis square
set(gcf,'unit','centimeter','Position',[7,4,13,5]);
box off