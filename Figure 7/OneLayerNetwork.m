%% One-population discrete attractor network - Connectivity, Tuning Curves

clc
clear
%% Define neurons

Nm = 300;
dthetam = 2*pi/Nm;
thetam = 0:dthetam:2*pi-dthetam;

Tau_syn = 0.01;
NEM = 1.5; thM = 0.1; sigM = 6.6; maxf = 100;
qm = @(x) maxf*(x-thM).^NEM./(sigM^NEM+(x-thM).^NEM).*(x>thM);
%% Generate memory circuit connectivity and background input

JE = 1;
JI = 0.17;
lambdaM = 0.2*pi;

alpha = 5e-4;
beta = 2.4e-3;

for i = 1:Nm
    thetaCurrent = thetam(i);
    JEModulated = JE*ConnModulation(thetaCurrent,alpha);
    JIModulated = JI*ConnModulation(thetaCurrent,beta);
    fE = JEModulated*exp(-((thetam-pi)/lambdaM).^2)/(2*pi);
    fI = JIModulated*exp(-((thetam-pi)/(0.6*pi)).^2)/(2*pi);
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

% For multiple input locations
NInputSample = 50;
dSample = 2*pi/NInputSample;
SampleInput = 0:dSample:2*pi;
NInputSample = NInputSample + 1;

%  % Compute preferred orietations
% 
%     TC = zeros(Nm,NInputSample);
%     PFTime = 500;
%     for jj = 1:NInputSample
%         mMemory_old = zeros(Nm,1);
%         MM_old = zeros(Nm,1);
%         I0 = (cos(thetam-pi)+1)/2;
%         I0 = circshift(I0,round(((-pi+SampleInput(jj))/dthetam)));
% 
%         for ii = 1:PFTime
%             MM_new = MM_old + 1/Tau_syn*dt*(-MM_old+mMemory_old);
%             MemoryInput = M*MM_old + IEc + I0'*(ii<(StimTime/dt));
%             mMemory_new = qm(MemoryInput);
%             mMemory_old = mMemory_new;
%             MM_old = MM_new;
%         end
%         TC(:,jj) = mMemory_new;
%     end
%     PF = ComputePF(TC,SampleInput);

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

WidthColor = '#4689A2';
FigOutDir = '';

f1 = figure;
figure(f1)
subplot(2,1,1)
plot(thetam,ME(1:20:end,:)/dthetam,'LineWidth',0.8,'Color','#1D2B79');
hold on
plot(thetam,MI(1:20:end,:)/dthetam,'LineWidth',0.7,'Color','#B81814');
hold off
box off
set(gca,'FontSize',12,'TickLength',[0.025,0.01],'TickDir','out','XTickLabelRotation',0);

yticks([0,0.1])
ylabel('Conn.','Interpreter','latex');
ylim([-0.05,0.18]);
xlim([0,2*pi]);
xticks(0:pi/2:2*pi);
xticklabels({'','','','',''});
subplot(2,1,2)
plot(thetam,M(1:20:end,:)/dthetam,'LineWidth',1);
xlim([0,2*pi]);
ylim([-0.02,0.18])
xticks(0:pi/2:2*pi);
yticks([0 0.1])
xticklabels({'0','','90','','180'});
set(gca,'FontSize',12,'TickLength',[0.025,0.01],'TickDir','out','XTickLabelRotation',0);
box off
xlabel('$\psi$ ($^\circ$)','Interpreter','latex');
set(gcf,'Unit','Centimeters','Position',[2,2,5,7.5]);

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
    ylim([0 60]);
    set(gca,'FontSize',10,'TickDir','out','TickLength',[0.04,.01],'LineWidth',.8, ...
        'XTickLabelRotation',0,'LooseInset',[0 0 0 0]);
end
set(gcf,'Unit','Centimeters','Position',[2,2,8,3]);

% f4 = figure;
% figure(f4)
% for i = 1:2
%     subplot(1,2,i)
%     PF = ComputePF(RateMemory(:,:,i),SampleInput);
%     histogram(PF,30,'FaceColor',[.5,.5,.5],'EdgeAlpha',0.5);
%     box off
% %     ylim([0 18]);
%     xlim([0 2*pi]);
%     yticks([0 15])
%     ylim([0,18])
%     if i == 1
%         ylabel('Counts');
%     end
%     set(gca,'FontSize',10,'TickDir','out','TickLength',[0.04,.01],'LineWidth',.8, ...
%         'XTickLabelRotation',0,'LooseInset',[0 0 0 0]);
%     xticks(0:pi/2:2*pi);
%     xticklabels({'0','','90','','180'});
%     xlabel('PF ($^\circ$)','Interpreter','latex');
% end
