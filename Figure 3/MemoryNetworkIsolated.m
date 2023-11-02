%% One-population discrete attractor network - Connectivity, Tuning Curves

clc
clear
%% Define neurons

Nm = 300;
dthetam = 2*pi/Nm; % We map 180 degs to 2*pi
thetam = 0:dthetam:2*pi-dthetam; % Neuron indices

Tau_syn = 0.01;
NEM = 1.5; thM = 0.1; sigM = 6.6; maxf = 100;
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
%% Euler integration

dt = 1e-3; % s
tmax = 4.5;
StimTime = 0.5;
PlotTime = StimTime + [0,4];
step = tmax/dt;

NInputSample = 50;
dSample = 2*pi/NInputSample;
SampleInput = 0:dSample:2*pi;
NInputSample = NInputSample + 1;

RateMemory = zeros(Nm,NInputSample,length(PlotTime));
tic
for jj = 1:NInputSample
    Iext = (cos(thetam-pi)+1)/2;
    I0(jj,:) = circshift(Iext,round(((-pi+SampleInput(jj))/dthetam))); %#ok<*SAGROW>
end
mMemory_old = zeros(Nm,NInputSample);
MM_old = zeros(Nm,NInputSample);
PlotTimeCounter = 1;

for ii = 1:step
    MM_new = MM_old + 1/Tau_syn*dt*(-MM_old+mMemory_old);
    MemoryInput = M*MM_old + IEc + I0'*(ii<(StimTime/dt));
    mMemory_new = qm(MemoryInput);
    if ii == round(PlotTime(PlotTimeCounter)/dt)
        RateMemory(:,:,PlotTimeCounter) = mMemory_new;
        if PlotTimeCounter < length(PlotTime)
            PlotTimeCounter = PlotTimeCounter + 1;
        end
    end
    mMemory_old = mMemory_new;
    MM_old = MM_new;
end
toc
%% Results

f1 = figure;
figure(f1)
plot(thetam,M(1:20:end,:)/dthetam,'LineWidth',1,'Color','#5C0B72');
xlim([0,2*pi]);
yticks([0 0.1])
ylim([-0.02 0.15])
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
set(gca,'FontSize',10,'TickLength',[0.025,0.01],'TickDir','out','XTickLabelRotation',0, ...
    'LooseInset',[0 0 0 0],'LineWidth',0.8);
box off
xlabel('$\psi$ ($^\circ$)','Interpreter','latex');
set(gcf,'Unit','Centimeters','Position',[2,2,5,3]);
ylabel('Conn.')

% Memory network tuning

f2 = figure;
figure(f2)
NeuronIdx = 1:20:Nm;
plot(SampleInput,RateMemory(NeuronIdx,:,2),'LineWidth',1);
set(gca,'FontSize',10);
xlim([0 2*pi]);
box off
ylim([0,40]);
xlabel('$\theta$ ($^\circ$)','Interpreter','latex');
set(gca,'FontSize',10,'TickDir','out','TickLength',[0.025,.01],'LineWidth',.8, ...
    'XTickLabelRotation',0,'LooseInset',[0 0 0 0]);
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
ylabel('FR (Hz)')
set(gcf,'Unit','Centimeters','Position',[2,2,5,3]);

f3 = figure;
figure(f3)
PF = ComputePF(RateMemory(:,:,2),SampleInput);
histogram(PF,30,'EdgeAlpha',0.5,'FaceColor',[0.5,0.5,0.5]);
box off
xlim([0 2*pi]);

ylabel('Counts','FontSize',10);
set(gca,'FontSize',10,'TickDir','out','TickLength',[0.025,.01],'LineWidth',.8, ...
    'XTickLabelRotation',0,'LooseInset',[0 0 0 0]);
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
xlabel('$\tilde{\theta}$','FontSize',10,'Interpreter','latex');
set(gcf,'Unit','Centimeters','Position',[2,2,5,2.5]);

WidthColor = '#5C0B72';
f4 = figure;
figure(f4)
Widths = zeros(1,Nm);
for k = 1:Nm
    Widths(k) = ComputeTuningWidth(SampleInput,RateMemory(k,:,2));
end
plot(thetam,Widths,'LineWidth',2,'Color',WidthColor);
set(gca,'FontSize',10,'TickDir','out','TickLength',[0.025,.01],'LineWidth',.8, ...
    'XTickLabelRotation',0,'LooseInset',[0 0 0 0]);
ylim([0 60]);
xlim([0 2*pi]);
box off
xlabel('\psi')
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
ylabel('Width ($^\circ$)','Interpreter','latex');
set(gcf,'Unit','Centimeters','Position',[2,2,5,2.5]);

%% Manifold

RSteadyState = RateMemory(:,:,2) - mean(RateMemory(:,:,2),2);
[U,S,V] = svd(RSteadyState');
fSVD2 = figure;

figure(fSVD2)
plot(V(:,1)'*RSteadyState,V(:,2)'*RSteadyState,'-','LineWidth',2,'Color', ...
    [0.7,0.7,0.7],'Marker','.','MarkerSize',14,'MarkerEdgeColor','#5C0B72');

axis square; axis off;
set(gca,'LooseInset',[0 0 0 0]);
set(gcf,'Unit','Centimeters','Position',[2,2,5,5]);

%%
fSVD3 = figure;
figure(fSVD3)

figure(fSVD3)
plot3(V(:,1)'*RSteadyState,V(:,2)'*RSteadyState,V(:,3)'*RSteadyState,'-', ...
    'LineWidth',3,'Color','#5C0B72');
axis square;
grid on
xlim([-220 220])
ylim([-220 220])
zlim([-220 220])
xticklabels({'','',''})
yticklabels({'','',''})
zticklabels({'','',''})
set(gcf,'Unit','Centimeters','Position',[2,2,5,5]);
view([-19 43])