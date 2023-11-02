%% Compare the full model and the simplified model

clear
close
SaveFlag = 0;
LoadFlag = 1;
%% Define parameters

if LoadFlag == 0
    seed = 123;
    rng(seed)

    SensoryNet.alpha = 0.03;
    SensoryNet.beta = 0;
    SensoryNet.Mode = 'EOnly';
    SensoryNet.N = 300;
    MemoryNet.N = 300;

    JForward = 0.1;
    nu = 0.17*pi;
    JBackward = 0.25;
    nufdbk = 0.17*pi;
    
    % Simulation parameters
    DynParams.NInputSample = round(SensoryNet.N/6)+1;
    DynParams.dSample = 2*pi/(DynParams.NInputSample-1);

    DynParams.dt = 1e-3;
    DynParams.Manifold_tmax = 1.5; % The manifold is parameterized at 1.5s after onset
    DynParams.StimTime = .5;  
    DynParams.NoiseTime = DynParams.Manifold_tmax;
    DynParams.RepTime = 3e3;  % The number of realizations
    DynParams.SimplifiedDecodeTime = [0.2,0.5,1];
    DynParams.FullModel_tmax = max(DynParams.SimplifiedDecodeTime) +...
                                    DynParams.Manifold_tmax;
    DynParams.FullDecodeTime = DynParams.SimplifiedDecodeTime + DynParams.Manifold_tmax;
    DynParams.PFTime = DynParams.Manifold_tmax;  % For estimation of PFs
    DynParams.Parallel = 16;
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
%% Evaluate v and D

    disp('Evaluating v_\theta and D_\theta...')
    tic
   [vtheta,Dtheta,ThetaLoc] = Compute_vtheta_Dtheta_new(SensoryNet,MemoryNet,DynParams);
    toc
    % vtheta, Dtheta, are in degrees
    % ThetaLoc is where the initial conditions we use for the simplified
    % model
    NFourier = 10;
    SampleInput = 0:DynParams.dSample:2*pi;
    [aD,bD] = FourierSeries(Dtheta*(2*pi/180)^2,NFourier,SampleInput); 
    [av,bv] = FourierSeries(vtheta*(2*pi/180),NFourier,SampleInput);  % here we convert 180 degrees to 2*pi
%% Dynamics of the simplified SDE

    disp('Solving the simplified SDE...');
    DecodedthetaSimplified = thetaSDEDynamics(av,bv,aD,bD,NFourier,DynParams,ThetaLoc);
%% Dynamics of the full SDE

    disp('Solving the full SDE...')
    [DecodedthetaFull,SensoryNet,MemoryNet,DynParams] = FullSDEDynamicsForComparison(SensoryNet,MemoryNet,DynParams);

    if SaveFlag
        save(['Comparisonalpha',strrep(num2str(SensoryNet.alpha),'.','p'),'beta', ...
            strrep(num2str(SensoryNet.beta),'.','p'),'PFTime',strrep(num2str(DynParams.PFTime),'.','p'),'Right']); %#ok<*UNRCH>
    end
    disp('----------Finished-----------');
end
%%
if LoadFlag
    clear
    close
    DataDir = '';
    load([DataDir,'Comparisonalpha0p05beta0PFTime1p5Right.mat']);
end

[BiasSimplified,BiasSimplifiedSEM,StdSimplified,StdSimplifiedSEM] = ...
        BiasStd(DecodedthetaSimplified,ThetaLoc);
[BiasFull,BiasFullSEM,StdFull,StdFullSEM] = BiasStd(DecodedthetaFull,ThetaLoc);


ColorSimplified = [0.2549,0.7137,0.7686;0.1137,0.5686,0.7529;0.1333,0.3686,0.6588];

ColorFull = [0.9961,0.6000,0.1608;0.9255,0.4392,0.0784;0.8000,0.2980,0.0078];
SampleInput = 0:DynParams.dSample:2*pi;

f1 = figure;
figure(f1)
subplot(1,2,1)
for Time = 1:size(BiasSimplified,1)
    Upper = BiasSimplified(Time,:) + BiasSimplifiedSEM(Time,:);
    Lower = BiasSimplified(Time,:) - BiasSimplifiedSEM(Time,:);
    fill([SampleInput, fliplr(SampleInput)],[Upper, fliplr(Lower)],ColorSimplified(Time,:),...
        'LineStyle','none','FaceAlpha',0.3);
    hold on
    plot(SampleInput,BiasSimplified(Time,:),'LineWidth',1.5,'Color',ColorSimplified(Time,:));
end

for Time = 1:size(BiasSimplified,1)
    Upper = BiasFull(Time,:) + BiasFullSEM(Time,:);
    Lower = BiasFull(Time,:) - BiasFullSEM(Time,:);
    fill([SampleInput, fliplr(SampleInput)],[Upper, fliplr(Lower)],ColorFull(Time,:),...
        'LineStyle','none','FaceAlpha',0.3);
    plot(SampleInput,BiasFull(Time,:),'LineWidth',1.5,'Color',ColorFull(Time,:));
end
hold off
% yline(0,'LineStyle','--','LineWidth',1.0);
set(gca,'FontSize',10,'LineWidth',.8,'TickLength',[0.025,0.01],'LooseInset',[0 0 0 0], ...
    'XTickLabelRotation',0);
xlabel('$\theta\; (^\circ)$','FontSize',12,'Interpreter','latex');
ylabel('Bias $(^\circ)$','FontSize',12,'Interpreter','latex');
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
xlim([0 2*pi]);
box off

subplot(1,2,2)
for Time = 1:size(BiasSimplified,1)
    Upper = StdSimplified(Time,:) + StdSimplifiedSEM(Time,:);
    Lower = StdSimplified(Time,:) - StdSimplifiedSEM(Time,:);
    fill([SampleInput, fliplr(SampleInput)],[Upper, fliplr(Lower)],ColorSimplified(Time,:),...
        'LineStyle','none','FaceAlpha',0.3);
    hold on
    plot(SampleInput,StdSimplified(Time,:),'LineWidth',1.5,'Color',ColorSimplified(Time,:));
end

for Time = 1:size(BiasSimplified,1)
    Upper = StdFull(Time,:) + StdFullSEM(Time,:);
    Lower = StdFull(Time,:) - StdFullSEM(Time,:);
    fill([SampleInput, fliplr(SampleInput)],[Upper, fliplr(Lower)],ColorFull(Time,:),...
        'LineStyle','none','FaceAlpha',0.3);
    plot(SampleInput,StdFull(Time,:),'LineWidth',1.5,'Color',ColorFull(Time,:));
end
hold off
set(gca,'FontSize',10,'LineWidth',.8,'TickLength',[0.025,0.01],'LooseInset',[0 0 0 0],'XTickLabelRotation',0);
xlabel('$\theta\; (^\circ)$','FontSize',12,'Interpreter','latex');
ylabel('SD $(^\circ)$','FontSize',12,'Interpreter','latex');
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
xlim([0 2*pi]);
ylim([0,1.5]);
box off
set(gcf,'unit','centimeter','Position',[2,2,9,4]);

%%
function [BiasMean,BiasSEM,Std,StdSEM] = BiasStd(Decodedtheta,ThetaLoc)
    NInputSample = size(Decodedtheta,2);
    RepTime = size(Decodedtheta,1);
    Bias = Decodedtheta-ThetaLoc;
    DecodeTimeLength = size(Decodedtheta,3);

    Bias(Bias>pi) = Bias(Bias>pi) - 2*pi;
    Bias(Bias<-pi) = Bias(Bias<-pi) + 2*pi;
    Bias = Bias/pi*180/2; % Note that 180 degrees were converted to 2*pi, so dividing by 2 is necessary
    BiasMean = zeros(DecodeTimeLength,NInputSample);
    BiasStdDv = BiasMean;
    BiasSEM = BiasStdDv;
    for Time = 1:DecodeTimeLength
        BiasMean(Time,:) = mean(Bias(:,:,Time),1);
        BiasStdDv(Time,:) = std(Bias(:,:,Time),0,1);
        BiasSEM(Time,:) = BiasStdDv(Time,:)/sqrt(RepTime);
    end
    Std = zeros(DecodeTimeLength,NInputSample);
    StdStdDv = Std;
    StdSEM = Std;
    for Time = 1:DecodeTimeLength
        Std(Time,:) = std(Bias(:,:,Time),0,1);
        MeanSquare = mean(Bias(:,:,Time),1).^2;
        StdSample = sqrt(Bias(:,:,Time).^2-MeanSquare);
        StdStdDv(Time,:) = std(StdSample,0,1);
        StdSEM(Time,:) = StdStdDv(Time,:)/sqrt(RepTime);
    end
end