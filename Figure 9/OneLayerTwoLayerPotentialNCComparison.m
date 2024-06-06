%% Compare the potential landscapes, tangent vector norms, and noise coefficients of the two-module and one-module model

clear
close
%% Define parameters for the two-layer network

seed = 123;
rng(seed)
SensoryNet.N = 300;
SensoryNet.alpha = 0.03;
SensoryNet.beta = 0.08;
SensoryNet.Mode = 'Both';

MemoryNet.N = 300;

JForward = 0.1;
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

%% Results of the two-layer network

tic
[vtheta,Dtheta,~] = Compute_vtheta_Dtheta_new(SensoryNet,MemoryNet,DynParams);
[dsdthetaNormTwoLayer,~] = ComputeTVNormTwoLayer(SensoryNet,MemoryNet,DynParams);
toc
SampleInput = 0:DynParams.dSample:2*pi;
TwoLayerPotential = EstimatePotential(vtheta,SampleInput);
TwoLayerNoiseCoef = sqrt(2*Dtheta);
%% Define the parameters of the memory single-layer network

MemoryNet.tau = 0.01;
MemoryNet.alpha = 5e-4;
MemoryNet.beta = 2.4e-3;
MemoryNet.N = 300;
MemoryNet.Mode = 'Both';
MemoryNet.NEM = 1.5; MemoryNet.thM = 0.1; MemoryNet.sigM = 6.6; MemoryNet.maxf = 100;
MemoryNet.q = @(x) MemoryNet.maxf*(x-MemoryNet.thM).^MemoryNet.NEM./( ...
    MemoryNet.sigM^MemoryNet.NEM+(x-MemoryNet.thM).^MemoryNet.NEM).*(x>MemoryNet.thM);

MemoryNet.JE = 1; MemoryNet.JI = 0.17;
MemoryNet.lambdaM = 0.2*pi;
MemoryNet.IEc = 0.6*ones(MemoryNet.N,1);
%% Results of the one-layer network

tic
MemoryNet = OneLayerRecurConn(MemoryNet);
[vtheta,Dtheta] = Compute_vtheta_Dtheta_OneLayer(MemoryNet,DynParams);
dsdthetaNormOneLayer = ComputeTVNormOneLayer(MemoryNet,DynParams);
toc
OneLayerPotential = EstimatePotential(vtheta,SampleInput);
OneLayerNoiseCoef = sqrt(2*Dtheta);
%% Plot results

TwoLayerColor = 'k';
OneLayerColor = '#c7111a';
FigOutDir = '';
f1 = figure;
figure(f1)
plot(SampleInput,TwoLayerPotential,'LineWidth',1,'Color',TwoLayerColor);
hold on
plot(SampleInput,OneLayerPotential,'LineWidth',1,'Color',OneLayerColor);
hold off
xlabel('$\theta~ (^\circ)$','Interpreter','latex');
xticks(0:pi/2:2*pi);
xlim([0 2*pi]);
ylim([0,max(OneLayerPotential)+1]);
xticklabels({'0','','90','','180'});
ylabel('Potential')
box off
set(gca,'FontSize',10,'LineWidth',0.8,'TickLength',[0.025,0.01],'LooseInset',[0 0 0 0],'TickDir','out');
set(gcf,'Units','Centimeters','Position',[2,2,5,4]);
% print(f1,'-depsc','-r300',[FigOutDir,'PotentialComparison.eps']);

f2 = figure;
figure(f2)
plot(SampleInput,TwoLayerNoiseCoef,'LineWidth',1,'Color','k');
hold on
plot(SampleInput,OneLayerNoiseCoef,'LineWidth',1,'Color',OneLayerColor);
hold off
xlabel('$\theta~ (^\circ)$','Interpreter','latex');
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
ylabel('Noise Coef.');
ylim([0,max(OneLayerNoiseCoef)+0.2]);
box off
set(gca,'FontSize',10,'LineWidth',0.8,'TickLength',[0.025,0.01],'LooseInset',[0 0 0 0],'TickDir','out');
set(gcf,'Units','Centimeters','Position',[2,2,5,4]);
% print(f2,'-depsc','-r300',[FigOutDir,'NoiseCoefComparison.eps']);

f3 = figure;
figure(f3)
plot(SampleInput,dsdthetaNormTwoLayer,'LineWidth',1,'Color','k');
hold on
plot(SampleInput,dsdthetaNormOneLayer,'LineWidth',1,'Color',OneLayerColor);
hold off
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
xlabel('$\theta~ (^\circ)$','Interpreter','latex');
ylabel('$\|ds/d\theta \| \; (\mathrm{deg}^{-1})$','Interpreter','latex');
ylim([0,max(dsdthetaNormTwoLayer)+2]);
box off
set(gca,'FontSize',10,'LineWidth',0.8,'TickLength',[0.025,0.01],'LooseInset',[0 0 0 0],'TickDir','out');
set(gcf,'Units','Centimeters','Position',[2,2,5,4]);
% print(f3,'-depsc','-r300',[FigOutDir,'dsdthetaComparison.eps']);


f4 = figure;
figure(f4)
plot(SampleInput,1./dsdthetaNormTwoLayer,'LineWidth',1,'Color','k');
hold on
plot(SampleInput,1./dsdthetaNormOneLayer,'LineWidth',1,'Color','#c30d23');
hold off
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
xlabel('$\theta~ (^\circ)$','Interpreter','latex');
ylabel('$\|ds/d\theta \|^{-1} \; (^\circ)$','Interpreter','latex');
ylim([0,max(1./dsdthetaNormOneLayer)*1.2]);
box off
set(gca,'FontSize',10,'LineWidth',0.8,'TickLength',[0.025,0.01],'LooseInset',[0 0 0 0],'TickDir','out');
set(gcf,'Units','Centimeters','Position',[2,2,5,4]);
% print(f4,'-depsc','-r300',[FigOutDir,'dsdthetaISComparison.eps']);