%% Plot the connectivity of the sensory circuit. Plot the tuning curves of the sensory-memory model.

clear
close
LoadFlag = 0;
SaveFlag = 0;
%% Define sensory neurons

if LoadFlag == 0
    Ns = 300;
    dthetas = 2*pi/Ns;
    Tau_syn = 0.01; % s
    thetas = 0:dthetas:2*pi-dthetas;
%% Generate sensory circuit connectivity

    JI = 0.35;
    JE = 0.6;
    alpha = 0.03;
    beta = 0.08;
    lambda = 0.36*pi;
    lambdaI = 1.1*pi;
    Conn = zeros(Ns,Ns);ConnI = Conn;ConnE = Conn;
    Mode = 'Both';

    for i = 1:Ns
        thetaCurrent = thetas(i);
        %   Modify gains
        if strcmp(Mode,'Both')
            JEModulated = JE*ConnModulation(thetaCurrent,alpha);
            JIModulated = JI*ConnModulation(thetaCurrent,beta);
            ConnI(i,:) = 1/2/pi*dthetas*circshift(-JIModulated*exp(-((thetas-pi)/lambdaI).^2),[0 round(-pi/dthetas-1+i)]);
        elseif strcmp(Mode,'EOnly') % Excitatory modulation
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

    %   % Input-output transfer function
    NEM = 1.5; thM = 0.1; sigM = 6.6;
    qm = @(x) maxf*(x-thM).^NEM./(sigM^NEM+(x-thM).^NEM).*(x>thM);
%% Generate memory circuit connectivity and background input

    JE = 1; JI = 0.17;
    lambdaM = 0.2*pi;

    ME = zeros(Nm,Nm);
    MI = ME;
    fE = JE*exp(-((thetam-pi)/lambdaM).^2)/(2*pi);
    fI = JI*exp(-((thetam-pi)/(pi*0.6)).^2)/(2*pi);
    for i = 1:Nm
        ME(i,:) = dthetam*circshift(fE,[0 -round(pi/dthetam)-1+i]);
        MI(i,:) = dthetam*circshift(fI,[0 -round(pi/dthetam)-1+i]);
    end

    M = ME-MI;

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
%% Define sensory input and Euler Integration

    C = 4;
    epsilon = 0.2;
    mu = 0.3*pi;

    dt = 1e-3; % ms

    % For multiple input locations
    NInputSample = 50;
    dSample = 2*pi/NInputSample;
    SampleInput = 0:dSample:2*pi;
    NInputSample = NInputSample + 1;
    StimTime = .5;
    % Compute preferred orietations

    PlotTime = StimTime+[0,4];
    tmax = PlotTime(end);
    step = tmax/dt;
    I0 = ExternalInput(NInputSample,Ns,C,epsilon,mu);

    RateSensory = zeros(Ns,NInputSample,length(PlotTime));
    RateMemory = zeros(Nm,NInputSample,length(PlotTime));
    % Synaptic variables
    SS_old = zeros(Ns,NInputSample);
    SM_old = zeros(Nm,NInputSample);
    MS_old = zeros(Ns,NInputSample);
    MM_old = zeros(Nm,NInputSample);
    % Firing rates
    mSensory_old = zeros(Ns,NInputSample);
    mMemory_old = zeros(Nm,NInputSample);
    PlotTimeCounter = 1;
    tic
    for ii = 1:step
        SS_new = SS_old+ 1/Tau_syn*dt*(-SS_old+mSensory_old);
        SM_new = SM_old + 1/Tau_syn*dt*(-SM_old+mMemory_old);
        MS_new = MS_old + 1/Tau_syn*dt*(-MS_old+mSensory_old);
        MM_new = MM_old + 1/Tau_syn*dt*(-MM_old+mMemory_old);
        SensoryInputI = Conn*SS_old + MBackward*SM_old+I0'*(ii<(StimTime/dt));

        mSensory_new = qs(SensoryInputI);
        MemoryInput = M*MM_old + IEc + MForward*MS_old;
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
%% Save data

    if SaveFlag
        SaveFlag = 0;
        disp('Saving the data.')
        save(['TuningCurveConnResultalpha',strrep(num2str(alpha),'.','p'), ...
            'beta',strrep(num2str(beta),'.','p'),'Tmax',num2str(tmax),'JBack', ...
            strrep(num2str(JBackward),'.','p'),'.mat']);
    end
end
%% Plot results
% Connectivity

if LoadFlag
    clear %#ok<*UNRCH> 
    load TuningCurveConnResultalpha0p03beta0p08Tmax4500JBack0p25.mat
end

WidthColor = '#4689A2';

if strcmp(Mode,'Both')
    f1 = figure;
    figure(f1)
    subplot(1,2,1)
    plot(thetas,ConnE(1:20:end,:)/dthetas,'LineWidth',1,'Color',[0.7,0.7,0.7]);
    hold on
    thetas_extended = [thetas(1:20:end),2*pi];
    EnvelopeE = max(ConnE(1:20:end,:)/dthetas,[],2);
    EnvelopeE = [EnvelopeE',EnvelopeE(1)];
    plot(thetas_extended,EnvelopeE,'LineWidth',2.5,'Color','#1D2B79');
    hold off
    box off
    set(gca,'FontSize',10,'TickLength',[0.025,0.01],'TickDir','out', ...
        'XTickLabelRotation',0,'LooseInset',[0 0 0 0],'LineWidth',0.8);

    yticks([0,0.1])
    ylabel('Conn.');
    ylim([-0,0.11]);
    xlim([0,2*pi]);
    xticks(0:pi/2:2*pi);
    xticklabels({'0','','90','','180'});
    subplot(1,2,2)
    plot(thetas,-ConnI(1:20:end,:)/dthetas,'LineWidth',1,'Color',[0.7 0.7 0.7]);
    hold on
    EnvelopeI = max(-ConnI(1:20:end,:)/dthetas,[],2);
    EnvelopeI = [EnvelopeI',EnvelopeI(1)];
    plot(thetas_extended,EnvelopeI,'LineWidth',2.5,'Color','#B81814');
    hold off
    xlim([0,2*pi]);
    xticks(0:pi/2:2*pi);
    ylim([0,0.065]);
    yticks([0,0.06])
    xticklabels({'0','','90','','180'});
    set(gca,'FontSize',10,'TickLength',[0.025,0.01],'TickDir','out', ...
        'XTickLabelRotation',0,'LooseInset',[0 0 0 0],'LineWidth',0.8);
    box off
    xlabel('$\psi$ ($^\circ$)','Interpreter','latex');
    set(gcf,'Unit','Centimeters','Position',[2,2,12,3]);

end

if strcmp(Mode,'EOnly')
    f1 = figure;
    figure(f1)
    subplot(1,2,1)
    plot(thetas,Conn(1:20:end,:)/dthetas,'LineWidth',1,'Color',[0.5,0.5,0.5]);
    Envelope = max(Conn(1:20:end,:)/dthetas,[],2);
    Envelope = [Envelope',Envelope(1)];
    hold on
    plot([thetas(1:20:end),2*pi],Envelope,'LineWidth',2,'Color','k');
    hold off
    xlim([0,2*pi]);
    yticks([-0.05 0 0.05])
    ylim([-0.06,0.05])
    xticks(0:pi/2:2*pi);
    xticklabels({'0','','90','','180'});
    set(gca,'FontSize',10,'TickLength',[0.025,0.01],'TickDir','out','XTickLabelRotation',0, ...
        'LooseInset',[0 0 0 0]);
    box off
    xlabel('$\psi$ ($^\circ$)','Interpreter','latex');

    subplot(1,2,2)
    plot(thetas,M(1:20:end,:)/dthetam,'LineWidth',1,'Color',[0.5,0.5,0.5]);
    xlim([0,2*pi]);
    ylim([-0.02,0.14])
    xticks(0:pi/2:2*pi);
    xticklabels({'0','','90','','180'});
    set(gca,'FontSize',10,'TickLength',[0.025,0.01],'TickDir','out','XTickLabelRotation',0, ...
        'LooseInset',[0 0 0 0]);
    box off
    xlabel('$\psi$ ($^\circ$)','Interpreter','latex');
    set(gcf,'Unit','Centimeters','Position',[2,2,10,3]);

end

% Memory network tuning

f3 = figure;
figure(f3)
NeuronIdx = 1:20:Nm;
for i = 1:2
    subplot(1,2,i)
    plot(SampleInput,RateMemory(NeuronIdx,:,i),'LineWidth',1);
    if i == 1
        ylabel('FR (Hz)');
    end
    xlim([0 2*pi]);
    box off
    xticks(0:pi/2:2*pi);
    xticklabels({'0','','90','','180'});
    TimeIntoDelay = (PlotTime(i)-StimTime);
    title(['Time = ', num2str(TimeIntoDelay),'s'],'FontSize',8);
    xlabel('$\theta$ ($^\circ$)','Interpreter','latex');
    ylim([0 70]);
    set(gca,'FontSize',10,'TickDir','out','TickLength',[0.04,.01],'LineWidth',.8, ...
        'XTickLabelRotation',0,'LooseInset',[0 0 0 0]);
end
set(gcf,'Unit','Centimeters','Position',[2,2,8,3]);

f4 = figure;
figure(f4)
for i = 1:2
    subplot(1,2,i)
    PF = ComputePF(RateMemory(:,:,i),SampleInput);
    histogram(PF,30,'FaceColor',[.5,.5,.5],'EdgeAlpha',0.5);
    box off
    xlim([0 2*pi]);
    yticks([0 15])
    ylim([0,19])
    if i == 1
        ylabel('Counts');
    end
    set(gca,'FontSize',10,'TickDir','out','TickLength',[0.04,.01],'LineWidth',.8, ...
        'XTickLabelRotation',0,'LooseInset',[0 0 0 0]);
    xticks(0:pi/2:2*pi);
    xticklabels({'0','','90','','180'});
    xlabel('PF ($^\circ$)','Interpreter','latex');
end
set(gcf,'Unit','Centimeters','Position',[2,2,8,3]);


f5 = figure;
figure(f5)
for i = 1:2
    subplot(1,2,i)
    TuningWidth = zeros(1,Nm);
    for j = 1:Nm
        TuningWidth(j) = ComputeTuningWidth(SampleInput,RateMemory(j,:,i));
    end
    plot(thetam,TuningWidth,'LineWidth',2,'Color',WidthColor);
    set(gca,'FontSize',10,'TickDir','out','TickLength',[0.04,.01],'LineWidth',.8, ...
        'XTickLabelRotation',0,'LooseInset',[0 0 0 0]);
    xlabel('$\psi$ ($^\circ$)','Interpreter','latex');
    if i == 1
        ylabel('Width ($^\circ$)','Interpreter','latex');
    end
    xlim([0,2*pi]);
    xticks(0:pi/2:2*pi);
    xticklabels({'0','','90','','180'});
    ylim([0,80]);
    yticks([0,80])
    box off
end
set(gcf,'Unit','Centimeters','Position',[2,2,8,3]);

% Sensory network tuning

f6 = figure;
figure(f6)
NeuronIdx = 1:20:Nm;
for i = 1:2
    subplot(1,2,i)
    plot(SampleInput,RateSensory(NeuronIdx,:,i),'LineWidth',1);
    if i == 1
        ylabel('FR (Hz)');
    end
    xlim([0 2*pi]);
    box off
    xticks(0:pi/2:2*pi);
    xticklabels({'0','','90','','180'});
    TimeIntoDelay = (PlotTime(i)-StimTime);
    title(['Time = ', num2str(TimeIntoDelay),'s'],'FontSize',8);
    xlabel('$\theta$ ($^\circ$)','Interpreter','latex');
    if i == 1
        ylim([0 40]);
    else
        ylim([0 6]);
    end
    set(gca,'FontSize',10,'TickDir','out','TickLength',[0.04,.01],'LineWidth',.8, ...
        'XTickLabelRotation',0,'LooseInset',[0 0 0 0]);
end
set(gcf,'Unit','Centimeters','Position',[2,2,8,3]);

f7 = figure;
figure(f7)
for i = 1:2
    subplot(1,2,i)
    PF = ComputePF(RateSensory(:,:,i),SampleInput);
    histogram(PF,30,'FaceColor',[.5,.5,.5],'EdgeAlpha',0.5);
    box off
    %     ylim([0 18]);
    xlim([0 2*pi]);
    yticks([0 15])
    ylim([0,19])
    if i == 1
        ylabel('Counts');
    end
    set(gca,'FontSize',10,'TickDir','out','TickLength',[0.04,.01],'LineWidth',.8, ...
        'XTickLabelRotation',0,'LooseInset',[0 0 0 0]);
    xticks(0:pi/2:2*pi);
    xticklabels({'0','','90','','180'});
    xlabel('PF ($^\circ$)','Interpreter','latex');
end
set(gcf,'Unit','Centimeters','Position',[2,2,8,3]);

%
f8 = figure;
figure(f8)
for i = 1:2
    subplot(1,2,i)
    TuningWidth = zeros(1,Nm);
    for j = 1:Nm
        TuningWidth(j) = ComputeTuningWidth(SampleInput,RateSensory(j,:,i));
    end
    plot(thetam,TuningWidth,'LineWidth',2,'Color',WidthColor);
    set(gca,'FontSize',10,'TickDir','out','TickLength',[0.04,.01],'LineWidth',.8, ...
        'XTickLabelRotation',0,'LooseInset',[0 0 0 0]);
    xlabel('$\psi$ ($^\circ$)','Interpreter','latex');
    if i == 1
        ylabel('Width ($^\circ$)','Interpreter','latex');
    end
    xlim([0,2*pi]);
    xticks(0:pi/2:2*pi);
    xticklabels({'0','','90','','180'});
    ylim([0,55]);
    yticks([0,50])
    box off
end
set(gcf,'Unit','Centimeters','Position',[2,2,8,3]);