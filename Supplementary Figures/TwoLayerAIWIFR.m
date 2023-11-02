%% Asymmetry index and width index -- associated with Figure 4

clear
close
LoadFlag = 0;
SaveFlag = 0;
%% Define sensory neurons

if LoadFlag == 0
    Ns = 300;
    Tau = 0.01; % s
    dthetas = 2*pi/Ns;
    thetas = 0:dthetas:2*pi-dthetas;
%% Generate sensory circuit connectivity

    JI = 0.35;
    JE = 0.6;
    alpha = 0.04;
    beta = 0;
    lambda = 0.36*pi;
    lambdaI = 1.1*pi;
    Conn = zeros(Ns,Ns);ConnI = Conn;ConnE = Conn;
    Mode = 'EOnly';

    TMSTime = 1000;
    for i = 1:Ns
        thetaCurrent = thetas(i);
        %   Modify gains
        if strcmp(Mode,'Both')
            JEModulated = JE*ConnModulation(thetaCurrent,alpha);
            JIModulated = JI*ConnModulation(thetaCurrent,beta);
            ConnI(i,:) = 1/2/pi*dthetas*circshift(-JIModulated*exp(-((thetas-pi)/lambdaI).^2),[0 round(-pi/dthetas-1+i)]);
        elseif strcmp(Mode,'EOnly')
            JEModulated = JE*ConnModulation(thetaCurrent,-alpha);
            JIModulated = JI*ConnModulation(thetaCurrent,0);
            ConnI(i,:) = 1/2/pi*dthetas*circshift(-JIModulated,[0 round(-pi/dthetas-1+i)]);
        end
        ConnE(i,:) = 1/2/pi*dthetas*circshift(JEModulated*exp(-((thetas-pi)/lambda).^2),[0,round(-pi/dthetas-1+i)]);
        Conn = ConnI+ConnE;
    end

    % Transfer function
    NE = 2; th = 0.1; sig = 6; maxf = 100;
    qs = @(x) maxf*(x-th).^NE./(sig^NE+(x-th).^NE).*(x>th);
%% Define memory neurons

    Nm = 300;
    dthetam = 2*pi/Nm;
    thetam = 0:dthetam:2*pi-dthetam;

    % Input-output transfer function
    NEM = 1.5; thM = 0.1; sigM = 6.6; maxf = 100;
    qm = @(x) maxf*(x-thM).^NEM./(sigM^NEM+(x-thM).^NEM).*(x>thM);
%% Generate memory circuit connectivity and background input

    JEM = 1; JIM = 0.17;
    lambdaM = 0.2*pi;
    ME = zeros(Nm,Nm);
    MI = zeros(Nm,Nm);

    for ii = 1:Nm
        thetaCurrent = thetam(ii);
        fE = JEM*exp(-((thetam-pi)/lambdaM).^2)/(2*pi);
        fI = JIM*exp(-((thetam-pi)/(pi*0.6)).^2)/(2*pi);
        ME(ii,:) = dthetam*circshift(fE,[0 round(-pi/dthetam-1+ii)]);
        MI(ii,:) = dthetam*circshift(fI,[0 -round(pi/(dthetam))-1+ii]);
    end

    M = ME - MI;

    IEc = 0.6*ones(Nm,1);
%% Generate feedforward connectivity

    JForward = .10;
    nu = 0.17*pi;
    dthetas = 2*pi/Ns;
    thetas = 0:dthetas:2*pi-dthetas;
    MForward = zeros(Nm,Ns);
    fForward = 1/(2*pi)*JForward*exp(-((pi-thetas)/nu).^2);
    for i= 1:Nm
        MForward(i,:) = dthetas*circshift(fForward,[0 ...
            -round(pi/dthetas)-1+i]);
    end
%% Generate feedback connectivity

    JBackward = .25;
    nufdbk = 0.17*pi;
    dthetam = 2*pi/Nm;
    thetam = 0:dthetam:2*pi-dthetam;
    MBackward = zeros(Ns,Nm);
    fBackward = 1/(2*pi)*JBackward*exp(-((pi-thetam)/nufdbk).^2);
    for i= 1:Ns
        MBackward(i,:) = dthetam*circshift(fBackward,[0 ...
            -round(pi/dthetam)-1+i]);
    end
%% Define sensory input

    C = 4;
    epsilon = 0.2;
    mu = 0.3*pi;

    NInputSample = 50;
    dSample = 2*pi/NInputSample;
    SampleInput = 0:dSample:2*pi;
    NInputSample = NInputSample + 1;
    I0 = ExternalInput(NInputSample,Ns,C,epsilon,mu);
%% Euler integration

    dt = 1e-3; % s
    tmax = 5;
    StimTime = .5;
    PlotTime = .05:.05:tmax;
    step = round(tmax/dt);

    RateSensory = zeros(Ns,NInputSample,length(PlotTime));
    RateMemory = zeros(Nm,NInputSample,length(PlotTime));

    % Compute preferred orietation for decoding
    PFTime = 5; % Saturation

    % Synaptic variables
    SS_old = zeros(Ns,NInputSample);
    SM_old = zeros(Nm,NInputSample);
    MS_old = zeros(Ns,NInputSample);
    MM_old = zeros(Nm,NInputSample);
    % Firing rates
    mSensory_old = zeros(Ns,NInputSample);
    mMemory_old = zeros(Nm,NInputSample);
    for ii = 1:round(PFTime/dt)
        SS_new = SS_old+ 1/Tau*dt*(-SS_old+mSensory_old);
        SM_new = SM_old + 1/Tau*dt*(-SM_old+mMemory_old);
        MS_new = MS_old + 1/Tau*dt*(-MS_old+mSensory_old);
        MM_new = MM_old + 1/Tau*dt*(-MM_old+mMemory_old);
        SensoryInputI = Conn*SS_old + MBackward*SM_old+I0'*(ii<(PFTime/dt));

        mSensory_new = qs(SensoryInputI);
        MemoryInput = M*MM_old + IEc + MForward*MS_old;
        mMemory_new = qm(MemoryInput);
        SS_old = SS_new;
        SM_old = SM_new;
        MS_old = MS_new;
        MM_old = MM_new;
        mSensory_old = mSensory_new;
        mMemory_old = mMemory_new;
    end
    PFSensory = ComputePF(mSensory_new,SampleInput);
    PFMemory = ComputePF(mMemory_new,SampleInput);

    tic
    mSensory_old = zeros(Ns,NInputSample);
    mMemory_old = zeros(Nm,NInputSample);
    SS_old = zeros(Ns,NInputSample);
    SM_old = zeros(Nm,NInputSample);
    MS_old = zeros(Ns,NInputSample);
    MM_old = zeros(Nm,NInputSample);
    PlotTimeCounter = 1;

    for ii = 1:step
        SS_new = SS_old+ 1/Tau*dt*(-SS_old+mSensory_old);
        SM_new = SM_old + 1/Tau*dt*(-SM_old+mMemory_old);
        MS_new = MS_old + 1/Tau*dt*(-MS_old+mSensory_old);
        MM_new = MM_old + 1/Tau*dt*(-MM_old+mMemory_old);
        SensoryInputI = Conn*SS_old + MBackward*SM_old+I0'*(ii<(StimTime/dt));
        mSensory_new = qs(SensoryInputI);
        MemoryInput = M*MM_old + IEc + MForward*MS_old*(ii<(TMSTime/dt));
        mMemory_new = qm(MemoryInput);
        if ii == round(PlotTime(PlotTimeCounter)/dt)
            RateSensory(:,:,PlotTimeCounter) = mSensory_new;
            RateMemory(:,:,PlotTimeCounter) = mMemory_new;
            if PlotTimeCounter < length(PlotTime)
                PlotTimeCounter = PlotTimeCounter + 1;
            end
        end

        SS_old = SS_new;
        SM_old = SM_new;
        MS_old = MS_new;
        MM_old = MM_new;
        mSensory_old = mSensory_new;
        mMemory_old = mMemory_new;
    end
    toc
%% Data analysis

    PlotNeuronInd = round(22.5/180*Ns); % Define sample sensory neuron 

    % Compute maximum firing rate of a sample sensory network neuron

    FRMaxSensory = zeros(1,length(PlotTime));
    for i = 1:length(PlotTime)
        FRMaxSensory(i) = max(RateSensory(PlotNeuronInd,:,i));
    end

    FRMaxMemory = zeros(1,length(PlotTime));
    for i = 1:length(PlotTime)
        FRMaxMemory(i) = max(RateMemory(PlotNeuronInd,:,i));
    end
% Compute asymmetry index

    AsyIndMemory = zeros(1,length(PlotTime));
    for i = 1:length(PlotTime)
        AsyIndMemory(i) = AsymmetryIndex(SampleInput,RateMemory(PlotNeuronInd,:,i));
    end

    AsyIndSensory = zeros(1,length(PlotTime));
    for i = 1:length(PlotTime)
        AsyIndSensory(i) = AsymmetryIndex(SampleInput,RateSensory(PlotNeuronInd,:,i));
    end
% Compute sensory biases

    BiasSensory = zeros(1,length(PlotTime));
    for Time = 1:length(PlotTime)
        for kk = 1:NInputSample
            BiasSensory(kk,Time) = PVDecoder(PFSensory,RateSensory(:,kk,Time))-SampleInput(kk);
        end
    end
    BiasSensory(BiasSensory>pi) = BiasSensory(BiasSensory>pi) - 2*pi;
    BiasSensory(BiasSensory<-pi) = BiasSensory(BiasSensory<-pi) + 2*pi;

    BiasSensory = 0.5*BiasSensory/pi*180;
    SampleBiasSensory = BiasSensory(round(NInputSample/8),:);
% Compute memory biases

    BiasMemory = zeros(1,length(PlotTime));
    for Time = 1:length(PlotTime)
        for k = 1:NInputSample
            BiasMemory(k,Time) = PVDecoder(PFMemory,RateMemory(:,k,Time))-SampleInput(k);
        end
    end
    BiasMemory(BiasMemory>pi) = BiasMemory(BiasMemory>pi) - 2*pi;
    BiasMemory(BiasMemory<-pi) = BiasMemory(BiasMemory<-pi) + 2*pi;

    BiasMemory = 0.5*BiasMemory/pi*180;

    SampleBiasMemory = BiasMemory(round(NInputSample/8),:);
% Compute tuning width index

    TuningWidthSensory = zeros(length(PlotTime),Ns);
    WidthContrastSensory = zeros(1,length(PlotTime));
    for i = 1:length(PlotTime)
        for j = 1:Ns
            TuningWidthSensory(i,j) = ComputeTuningWidth(SampleInput,RateSensory(j,:,i));
        end
        MaxWidth = max(TuningWidthSensory(i,:));
        MinWidth = min(TuningWidthSensory(i,:));
        WidthContrastSensory(i) = (MaxWidth-MinWidth)/(MaxWidth+MinWidth);
    end

    TuningWidthMemory = zeros(length(PlotTime),Nm);
    WidthContrastMemory = zeros(1,length(PlotTime));
    for i = 1:length(PlotTime)
        for j = 1:Nm
            TuningWidthMemory(i,j) = ComputeTuningWidth(SampleInput,RateMemory(j,:,i));
        end
        MaxWidth = max(TuningWidthMemory(i,:));
        MinWidth = min(TuningWidthMemory(i,:));
        WidthContrastMemory(i) = (MaxWidth-MinWidth)/(MaxWidth+MinWidth);
    end
%% Save data

    if SaveFlag
        SaveFlag = 0;
        disp('Saving the data.')
        DataDir = '';
        save([DataDir,'AsymmetryTuningWidthFR','Mode',Mode,'alpha',strrep(num2str(alpha),'.','p'),'beta',strrep(num2str(beta),'.','p')]);
    end
end
% Plot bias increase

if LoadFlag
    close
    clear
     DataDir = '';
    load([DataDir,'AsymmetryTuningWidthFRModeBothalpha0p03beta0p08.mat']);
end

MemColor = '#5C0B72';
SenColor = '#448983';

f1 = figure;
figure(f1)
TimeInSecond = PlotTime;
plot(TimeInSecond,SampleBiasSensory,'LineWidth',2,'Color',SenColor);
hold on
plot(TimeInSecond,SampleBiasMemory,'LineWidth',2,'Color',MemColor);
hold off
xline(StimTime,'Color',[0.3,0.3,0.3],'LineWidth',1.2);
ylim([0,max(SampleBiasMemory)+0.05]);
ylabel('Bias (deg)','FontSize',10);
box off
set(gca,'FontSize',10,'LooseInset',[0 0 0 0],'TickLength',[0.025,0.01],'TickDir','out','LineWidth',0.8);
xlim([0,max(TimeInSecond)]);
xline(TMSTime,'Color',[0.3,0.3,0.3],'LineWidth',1.2);
xlabel('Time (s)','FontSize',10);
set(gcf,'Units','Centimeters','Position',[2,2,4,3]);

% Plot asymmetry index

f2 = figure;
figure(f2)
plot(TimeInSecond,AsyIndSensory,'LineWidth',2,'Color',SenColor);
hold on
plot(TimeInSecond,AsyIndMemory,'LineWidth',2,'Color',MemColor);
ylabel('AI (a.u.)','FontSize',10);
xline(StimTime,'Color',[0.3,0.3,0.3],'LineWidth',1.2);
box off
hold off
set(gca,'FontSize',10,'LooseInset',[0 0 0 0],'TickLength',[0.025,0.01],'TickDir','out','LineWidth',0.8);
xlabel('Time (s)','FontSize',10);
xlim([0,max(TimeInSecond)]);
xline(TMSTime,'Color',[0.3,0.3,0.3],'LineWidth',1.2);
ylim([0,1]);
yticks([0,1]);
set(gcf,'Units','Centimeters','Position',[2,2,4,3]);

% Plot firing rate change

f3 = figure;
figure(f3)
TimeInSecond = PlotTime;
plot(TimeInSecond,FRMaxSensory,'LineWidth',2,'Color',SenColor);
hold on
plot(TimeInSecond,FRMaxMemory,'LineWidth',2,'Color',MemColor);
hold off
xline(StimTime,'Color',[0.3,0.3,0.3],'LineWidth',1.2);
ylim([0,max(FRMaxSensory)+2]);
ylabel('FR (Hz)','FontSize',10);
box off
set(gca,'FontSize',10,'LooseInset',[0 0 0 0],'TickLength',[0.025,0.01],'TickDir','out','LineWidth',0.8);
xlim([0,max(TimeInSecond)]);
xline(TMSTime,'Color',[0.3,0.3,0.3],'LineWidth',1.2);
xlabel('Time (s)','FontSize',10);
set(gcf,'Units','Centimeters','Position',[2,2,4,3]);

% Plot tuning widths change

f6 = figure;
figure(f6)
plot(TimeInSecond,WidthContrastSensory,'LineWidth',2,'Color',SenColor);
hold on
plot(TimeInSecond,WidthContrastMemory,'LineWidth',2,'Color',MemColor);
hold off
xlim([0,max(TimeInSecond)]);
xline(StimTime,'Color',[0.3,0.3,0.3],'LineWidth',1.2);
xlabel('Time (s)','FontSize',10);
xline(TMSTime,'Color',[0.3,0.3,0.3],'LineWidth',1.2);
set(gca,'FontSize',10,'LooseInset',[0 0 0 0],'TickLength',[0.025,0.01],'TickDir','out','LineWidth',0.8);
box off
ylabel('WI (a.u.)','FontSize',10);
set(gcf,'Units','Centimeters','Position',[2,2,4,3]);