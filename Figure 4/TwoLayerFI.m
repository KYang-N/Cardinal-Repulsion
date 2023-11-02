%% Simulation of the two-layer network with Poisson noise
%% Fisher information

clear
close
SaveFlag = 0;
LoadFlag = 1;
%% Define parameters

if LoadFlag == 0
    seed = 233;
    rng(seed)

    SensoryNet.N = 300;
    SensoryNet.alpha = 0.04;
    SensoryNet.beta = 0;
    SensoryNet.Mode = 'EOnly';

    MemoryNet.N = 300;

    JForward = .10;
    nu = 0.17*pi;
    JBackward = .25;
    nufdbk = 0.17*pi;

    % Simulation parameters
    DynParams.AddNoise = 1;
    DynParams.dt = 1e-3;
    DynParams.StimTime = 0.5;
    DynParams.RepTime = 3e3;
    DynParams.FullDecodeTime = [1,2.5,4] + DynParams.StimTime;
    DynParams.Parallel = 4;  % Number of parallel workers
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
    [FI,SensoryNet,MemoryNet,DynParams] = TwoLayerNetworkFI(SensoryNet,MemoryNet,DynParams);
    delete(p);
    if SaveFlag
        disp('Saving the data.')
        save(['TwoLayerFIResultsalpha',strrep(num2str(SensoryNet.alpha),'.','p'), ...
            'beta',strrep(num2str(SensoryNet.beta),'.','p'),'Mode',SensoryNet.Mode,'.mat']);
    end
end
%%
if LoadFlag
    clear %#ok<*UNRCH>
    DataDir = '';
    load([DataDir,'TwoLayerFIResultsalpha0p04beta0ModeEonly.mat']);
end
close

Color = [0.6667 0.6667 0.6667;0.3333 0.3333 0.3333;0 0 0];
%% Plot FI

SampleInput = 0:DynParams.dSample:2*pi;
FigOutDir = '';
f = figure;
figure(f)
for Time = 1:3
    %     plot(SampleInput,FI(Time,:)/max(FI(1,:)),'Color',Color(Time,:),'LineWidth',1.5);
    plot(SampleInput,FI(Time,:),'Color',Color(Time,:),'LineWidth',1.5);
    hold on
end
hold off
xlabel('$\theta$ ($^\circ$)','Interpreter','latex');
ylabel('FI (deg^{-2})');
set(gca,'FontSize',10,'LineWidth',1,'TickLength',[0.025,0.01],'TickDir','out','XTickLabelRotation',0, ...
    'LooseInset',[0 0 0 0]);
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
box off
set(gcf,'Unit','Centimeters','Position',[2,2,4,4]);
