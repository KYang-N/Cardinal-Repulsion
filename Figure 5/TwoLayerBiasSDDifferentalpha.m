%% Simulation of the two-layer network with Poisson noise
%% Bias, SD -- Different \alpha

clear
close
SaveFlag = 0;
LoadFlag = 1;
%% Define parameters

if LoadFlag == 0
    seed = 123;
    rng(seed)

    SensoryNet.N = 300;
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
    DynParams.RepTime = 1e3;
    DynParams.PFTime = 5;
    DynParams.NInputSample = 51; % The first input orientation is 0, the last is 2*pi
    DynParams.FullDecodeTime = 4 + DynParams.StimTime;
    DynParams.Decoder = 'Aware';
    DynParams.Parallel = 10;
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

    alphaRange = [0.03 0.04 0.05];
    DecodedValue = zeros(DynParams.RepTime,DynParams.NInputSample,length(alphaRange));
    for ii = 1:length(alphaRange)
        SensoryNet.alpha = alphaRange(ii);
        [DecodedthetaFull,SesoryNet,MemoryNet,DynParams] = FullSDEDynamics(SensoryNet,MemoryNet,DynParams);
        DecodedValue(:,:,ii) = reshape(DecodedthetaFull,DynParams.RepTime,DynParams.NInputSample);
    end
    if SaveFlag
        disp('Saving the data.')
        save(['PoissonNoiseResultsDifferentalpha','Mode',SensoryNet.Mode,'.mat']);
    end
end
%%
if LoadFlag
    clear %#ok<*UNRCH> 
    DataDir = '';
    load([DataDir,'PoissonNoiseResultsDifferentalphaModeEOnly.mat']);
end
close
Color = [0.6667 0.6667 0.6667;0.3333 0.3333 0.3333;0 0 0];
%% Plot bias, SD, and discriminability pattern

SampleInput = 0:DynParams.dSample:2*pi;
FigOutDir = '';
NInputSample = size(DecodedValue,2);
RepTime = DynParams.RepTime;

Bias = DecodedValue-SampleInput;

Bias(Bias>pi) = Bias(Bias>pi) - 2*pi;
Bias(Bias<-pi) = Bias(Bias<-pi) + 2*pi;
Bias = Bias/pi*180/2; % In the simulation, 180 degrees are converted to 2*pi

BiasMean = zeros(length(alphaRange),NInputSample);
StdDv = BiasMean;
SEM = StdDv;

for Time = 1:length(alphaRange)
    BiasMean(Time,:) = mean(Bias(:,:,Time),1);
    StdDv(Time,:) = std(Bias(:,:,Time),0,1);
    SEM(Time,:) = StdDv(Time,:)/sqrt(RepTime);
end

f1 = figure;
figure(f1)
subplot(1,2,1)
for Time = 1:length(alphaRange)
    Upper = BiasMean(Time,:) + SEM(Time,:);
    Lower = BiasMean(Time,:) - SEM(Time,:);
    fill([SampleInput, fliplr(SampleInput)],[Upper, fliplr(Lower)],Color(Time,:),...
        'LineStyle','none','FaceAlpha',0.3);
    hold on
    plot(SampleInput,BiasMean(Time,:),'LineWidth',1.5,'Color',Color(Time,:));
end
hold off
yline(0,'LineStyle','--','Color','k','LineWidth',1.0);
xlabel('\theta (deg)');
ylabel('Bias (deg)');
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
xlim([0 2*pi]);
ylim([min(BiasMean(end,:)-1),max(BiasMean(end,:)+1)])
set(gca,'FontSize',10,'LineWidth',0.8,'TickLength',[0.025,0.01],'TickDir','out','XTickLabelRotation',0, ...
    'LooseInset',[0 0 0 0]);
box off


%% Std
Std = zeros(length(alphaRange),NInputSample);
for Time = 1:length(alphaRange)
    Std(Time,:) = std(Bias(:,:,Time),0,1);
    MeanSquare = mean(Bias(:,:,Time),1).^2;
    StdSample = sqrt(Bias(:,:,Time).^2-MeanSquare);
    StdDv(Time,:) = std(StdSample,0,1);
    SEM(Time,:) = StdDv(Time,:)/sqrt(RepTime);
end

figure(f1)
subplot(1,2,2)
for Time = 1:length(alphaRange)
    Upper = Std(Time,:) + SEM(Time,:);
    Lower = Std(Time,:) - SEM(Time,:);
    fill([SampleInput, fliplr(SampleInput)],[Upper, fliplr(Lower)],Color(Time,:),...
        'LineStyle','none','FaceAlpha',0.3);
    hold on
    plot(SampleInput,Std(Time,:),'LineWidth',1.5,'Color',Color(Time,:));
end
hold off
set(gca,'FontSize',10,'LineWidth',0.8,'TickLength',[0.025,0.01],'TickDir','out','XTickLabelRotation',0, ...
    'LooseInset',[0 0 0 0]);

ylim([1.8 max(Std(end,:))+0.1]);
xlim([0 2*pi]);
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
box off
xlabel('\theta (deg)');
ylabel('SD (deg)');
set(f1,'Units','Centimeters','Position',[2,2,10,5]);