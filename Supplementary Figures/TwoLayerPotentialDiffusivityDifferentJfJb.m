%% Using projection to estimate the drift velocity and noise coefficient of \theta for different Jf and Jb

clear
close
%% Define parameters

seed = 123;
rng(seed)
SensoryNet.N = 300;
% SensoryNet.beta = 0.09;
SensoryNet.alpha = 0.04;
SensoryNet.Mode = 'EOnly';

MemoryNet.N = 300;

JForwardRange = [0.07,0.1,0.13];
nu = 0.17*pi;
JBackward = 0.25;
nufdbk = 0.17*pi;

% Simulation parameters
DynParams.NInputSample = round(SensoryNet.N/6)+1;
DynParams.dSample = 2*pi/(DynParams.NInputSample-1);
DynParams.AddNoise = 0;
DynParams.dt = 1e-3;
DynParams.StimTime = 0.5;
DynParams.Manifold_tmax = 1.5;
DynParams.PFTime = 1.5;
%% Generate feedback connectivity

dthetam = 2*pi/MemoryNet.N;
thetam = 0:dthetam:2*pi-dthetam;
MemoryNet.MBackward = zeros(SensoryNet.N,MemoryNet.N);
fBackward = 1/(2*pi)*JBackward*exp(-((pi-thetam)/nufdbk).^2);
for i= 1:SensoryNet.N
    MemoryNet.MBackward(i,:) = dthetam*circshift(fBackward,[0 ...
        -round(pi/dthetam)-1+i]);
end
%% Different Jf

vtheta = zeros(length(JForwardRange),DynParams.NInputSample);
Dtheta = zeros(length(JForwardRange),DynParams.NInputSample);
tic
for ii = 1:length(JForwardRange)
    JForward = JForwardRange(ii);
    % Generate feedforward connectivity
    dthetas = 2*pi/SensoryNet.N;
    thetas = 0:dthetas:2*pi-dthetas;
    SensoryNet.MForward = zeros(MemoryNet.N,SensoryNet.N);
    fForward = 1/(2*pi)*JForward*exp(-((pi-thetas)/nu).^2);
    for i= 1:MemoryNet.N
        SensoryNet.MForward(i,:) = dthetas*circshift(fForward,[0 ...
            -round(pi/dthetas)-1+i]);
    end
    [vtheta(ii,:),Dtheta(ii,:),~] = Compute_vtheta_Dtheta_new(SensoryNet,MemoryNet,DynParams);
end
toc
%%
SampleInput = 0:DynParams.dSample:2*pi;
Color = [0.6667 0.6667 0.6667;0.3333 0.3333 0.3333;0 0 0];

f1 = figure;
figure(f1)
for ii = 1:length(JForwardRange)
    subplot(1,3,1)
    Potential = EstimatePotential(vtheta(ii,:),SampleInput);
    plot(SampleInput,Potential,'LineWidth',2,'Color',Color(ii,:));
    hold on

    subplot(1,3,2)
    plot(SampleInput,vtheta(ii,:),'LineWidth',2,'Color',Color(ii,:));
    hold on

    subplot(1,3,3)
    plot(SampleInput,sqrt(2*Dtheta(ii,:)),'LineWidth',2,'Color',Color(ii,:));
    hold on
end
subplot(1,3,1)
xlabel('$\theta$ ($^\circ$)','Interpreter','latex');
ylabel('Potential (degrees^2/s)');
ylim([0 25]);
set(gca,'FontSize',10,'TickLength',[0.025,0.01],'XTickLabelRotation',0,'LineWidth',0.8,'TickDir','out');
xlim([0 2*pi]);
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
box off

subplot(1,3,2)
set(gca,'FontSize',10,'TickLength',[0.025,0.01],'XTickLabelRotation',0,'LineWidth',0.8,'TickDir','out');
box off
yline(0,'--','LineWidth',1.0);
xlabel('$\theta$ ($^\circ$)','Interpreter','latex');
ylabel('Velocity ($^\circ/\mathrm{s}$)','Interpreter','latex');
xticks(0:pi/2:2*pi);
ylim([-0.9,0.9]);
xticklabels({'0','','90','','180'});


subplot(1,3,3)
set(gca,'FontSize',10,'TickLength',[0.025,0.01],'XTickLabelRotation',0,'LineWidth',0.8,'TickDir','out');
xlim([0 2*pi]);
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
xlabel('$\theta$ ($^\circ$)','Interpreter','latex');
ylabel('Noise coef. ($^\circ$)','Interpreter','latex');
xlim([0 2*pi]);
ylim([0,1.4]);
box off
set(gcf,'Unit','Centimeters','Position',[2,2,13.5,4]);

%% Different Jb

JBackwardRange = [0.22,0.25,0.28];
JForward = 0.1;
% Generate feedforward connectivity
dthetas = 2*pi/SensoryNet.N;
thetas = 0:dthetas:2*pi-dthetas;
SensoryNet.MForward = zeros(MemoryNet.N,SensoryNet.N);
fForward = 1/(2*pi)*JForward*exp(-((pi-thetas)/nu).^2);
for i= 1:MemoryNet.N
    SensoryNet.MForward(i,:) = dthetas*circshift(fForward,[0 ...
        -round(pi/dthetas)-1+i]);
end

vtheta = zeros(length(JForwardRange),DynParams.NInputSample);
Dtheta = zeros(length(JForwardRange),DynParams.NInputSample);
tic
for ii = 1:length(JBackwardRange)
    JBackward= JBackwardRange(ii);
    % Generate feedback connectivity
    dthetam = 2*pi/MemoryNet.N;
    thetam = 0:dthetam:2*pi-dthetam;
    MemoryNet.MBackward = zeros(SensoryNet.N,MemoryNet.N);
    fBackward = 1/(2*pi)*JBackward*exp(-((pi-thetam)/nufdbk).^2);
    for i= 1:MemoryNet.N
        MemoryNet.MBackward(i,:) = dthetam*circshift(fBackward,[0 ...
            -round(pi/dthetam)-1+i]);
    end

    [vtheta(ii,:),Dtheta(ii,:),~] = Compute_vtheta_Dtheta_new(SensoryNet,MemoryNet,DynParams);
end
toc
%%
f2 = figure;
figure(f2)
for ii = 1:length(JBackwardRange)
    subplot(1,3,1)
    Potential = EstimatePotential(vtheta(ii,:),SampleInput);
    plot(SampleInput,Potential,'LineWidth',2,'Color',Color(ii,:));
    hold on

    subplot(1,3,2)
    plot(SampleInput,vtheta(ii,:),'LineWidth',2,'Color',Color(ii,:));
    hold on

    subplot(1,3,3)
    plot(SampleInput,sqrt(2*Dtheta(ii,:)),'LineWidth',2,'Color',Color(ii,:));
    hold on
end
subplot(1,3,1)
xlabel('$\theta$ ($^\circ$)','Interpreter','latex');
ylabel('Potential (degrees^2/s)');
ylim([0 30]);
set(gca,'FontSize',10,'TickLength',[0.025,0.01],'XTickLabelRotation',0,'LineWidth',0.8,'TickDir','out');
xlim([0 2*pi]);
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
box off

subplot(1,3,2)
set(gca,'FontSize',10,'TickLength',[0.025,0.01],'XTickLabelRotation',0,'LineWidth',0.8,'TickDir','out');
box off
yline(0,'--','LineWidth',1.0);
xlabel('$\theta$ ($^\circ$)','Interpreter','latex');
ylabel('Velocity ($^\circ/\mathrm{s}$)','Interpreter','latex');
xticks(0:pi/2:2*pi);
ylim([-1.2,1.2]);
xticklabels({'0','','90','','180'});


subplot(1,3,3)
set(gca,'FontSize',10,'TickLength',[0.025,0.01],'XTickLabelRotation',0,'LineWidth',0.8,'TickDir','out');
xlim([0 2*pi]);
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
xlabel('$\theta$ ($^\circ$)','Interpreter','latex');
ylabel('Noise coef. ($^\circ$)','Interpreter','latex');
xlim([0 2*pi]);
ylim([0,1.4]);
box off
set(gcf,'Unit','Centimeters','Position',[2,2,13.5,4]);