%% Simulation of the two-layer network and generate schematic of the potential well

clear
LoadFlag = 0;
SaveFlag = 0;
%% Define parameters

Ns = 300;
Tau = 0.01; % s
dthetas = 2*pi/Ns;
thetas = 0:dthetas:2*pi-dthetas;
%% Generate sensory circuit connectivity

JI = 0.35;
JE = 0.6;
alpha = 0.05;
beta = 0;
lambda = 0.36*pi;
lambdaI = 1.1*pi;
Conn = zeros(Ns,Ns);ConnI = Conn;ConnE = Conn;
Mode = 'EOnly';

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

NInputSample = 100;
dSample = 2*pi/NInputSample;
SampleInput = 0:dSample:2*pi;
NInputSample = NInputSample + 1;
I0 = ExternalInput(NInputSample,Ns,C,epsilon,mu);
%% Euler integration

dt = 1e-3; % s
tmax = 4;
StimTime = .5;
step = round(tmax/dt);

tic
mSensory_old = zeros(Ns,NInputSample);
mMemory_old = zeros(Nm,NInputSample);
SS_old = zeros(Ns,NInputSample);
SM_old = zeros(Nm,NInputSample);
MS_old = zeros(Ns,NInputSample);
MM_old = zeros(Nm,NInputSample);

for ii = 1:step
    SS_new = SS_old+ 1/Tau*dt*(-SS_old+mSensory_old);
    SM_new = SM_old + 1/Tau*dt*(-SM_old+mMemory_old);
    MS_new = MS_old + 1/Tau*dt*(-MS_old+mSensory_old);
    MM_new = MM_old + 1/Tau*dt*(-MM_old+mMemory_old);
    SensoryInputI = Conn*SS_old + MBackward*SM_old+I0'*(ii<(StimTime/dt));
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
toc
%% Memory manifold PCA

MemoryState = [SS_old;SM_old;MS_old;MM_old];
StateMean = mean(MemoryState,2);
MemoryState = MemoryState - StateMean;

%% Generate a few trajectories

if LoadFlag == 0
    dt = 1e-3; % s
    tmax = 20;
    StimTime = .5;
    step = round(tmax/dt);
    InputOrientation = (pi+3/2*pi)/2;
    rng(123);
    PlotTime = 0.5:0.1:tmax;
    Trajectory = zeros(4*Ns,length(PlotTime),5);
    MemoryFiringRate = zeros(Nm,length(PlotTime),5);
    I_ext = C*(1-2*epsilon+2*epsilon*exp(-((thetas-pi)/mu).^2));
    I0 = circshift(I_ext,round(((-pi+InputOrientation)/dthetas)));
    tic
    for ll = 1:5
        mSensory_old = zeros(Ns,1);
        mMemory_old = zeros(Nm,1);
        SS_old = zeros(Ns,1);
        SM_old = zeros(Nm,1);
        MS_old = zeros(Ns,1);
        MM_old = zeros(Nm,1);
        PlotTimeCounter = 1;

        for ii = 1:step

            SS_new = SS_old+ 1/Tau*dt*(-SS_old+mSensory_old)+ ...
                1/Tau*sqrt(dt)*sqrt(mSensory_old*dt).*randn(Ns,1);
            SM_new = SM_old + 1/Tau*dt*(-SM_old+mMemory_old)+...
                1/Tau*sqrt(dt)*sqrt(mMemory_old*dt).*randn(Nm,1);
            MS_new = MS_old + 1/Tau*dt*(-MS_old+mSensory_old)+...
                1/Tau*sqrt(dt)*sqrt(mSensory_old*dt).*randn(Ns,1);
            MM_new = MM_old + 1/Tau*dt*(-MM_old+mMemory_old)+...
                1/Tau*sqrt(dt)*sqrt(mMemory_old*dt).*randn(Nm,1);
            SensoryInputI = Conn*SS_old + MBackward*SM_old+I0'*(ii<(StimTime/dt));
            mSensory_new = qs(SensoryInputI);
            MemoryInput = M*MM_old + IEc + MForward*MS_old;
            mMemory_new = qm(MemoryInput);


            if ii == round(PlotTime(PlotTimeCounter)/dt)
                Trajectory(:,PlotTimeCounter,ll) = [SS_new;SM_new;MS_new;MM_new];
                MemoryFiringRate(:,PlotTimeCounter,ll) = mMemory_new;
                PlotTimeCounter = PlotTimeCounter + 1;
            end

            SS_old = SS_new;
            SM_old = SM_new;
            MS_old = MS_new;
            MM_old = MM_new;
            mSensory_old = mSensory_new;
            mMemory_old = mMemory_new;
        end
    end
    toc
    if SaveFlag
        save(['Trajectoryalpha',strrep(num2str(alpha),'.','p'),'.mat'],'Trajectory','MemoryFiringRate');
    end
end
if LoadFlag
    DataDir = '';
    load([DataDir,'Trajectoryalpha0p05.mat']);
end

[U,S,V] = svd(MemoryState');

ManifoldColor = '#0087A9';

f1 = figure;
figure(f1)
PCaxis1 = V(:,1);
PCaxis2 = V(:,2);
PC1 = V(:,1)'*MemoryState;
PC2 = V(:,2)'*MemoryState;
MeanPC1 = PCaxis1'*StateMean;
MeanPC2 = PCaxis2'*StateMean;
plot(PC1+MeanPC1,PC2+MeanPC2,'LineWidth',4,'Color',ManifoldColor);
box off
set(gca,'FontSize',10,'LineWidth',.8,'LooseInset',[0 0 0 0],'TickLength',[0.025,0.01],'TickDir','out');
hold on
axis off
% Project vector field
Step = 120;
Pad = 20;
[u,v] = meshgrid(-450-MeanPC1-Pad:Step:500+PC1+Pad,-450-MeanPC2-Pad:Step:500+MeanPC2+Pad);
VectorFieldPC1 = u;
VectorFieldPC2 = v;
for ii = 1:size(u,1)
    for jj = 1:size(u,2)
        State = u(ii,jj)*PCaxis1 + v(ii,jj)*PCaxis2+StateMean;
        SS = State(1:Ns);
        SM = State(Ns+1:2*Ns);
        MS = State(2*Ns+1:3*Ns);
        MM = State(3*Ns+1:4*Ns);
        SensoryInput = Conn*SS + MBackward*SM;
        MemoryInput = M*MM + IEc + MForward*MS;
        mSensory = qs(SensoryInput);
        mMemory = qm(MemoryInput);
        VF = [-SS+mSensory;-SM+mMemory;-MS+mSensory;-MM+mMemory];
        VectorFieldPC1(ii,jj) = PCaxis1'*VF;
        VectorFieldPC2(ii,jj) = PCaxis2'*VF;
    end
end
xlim([-450,450]);
ylim([-450,450]);
quiver(u,v,VectorFieldPC1,VectorFieldPC2,'Color','k','LineWidth',0.8,'AutoScaleFactor',1);
hold off
axis square
set(gcf,'Units','Centimeters','Position',[2,2,6,6]);
hold on
for ll = 1:5
    TrajectoryPC1 = PCaxis1'*(Trajectory(:,:,ll)-StateMean)+MeanPC1;
    TrajectoryPC2 = PCaxis2'*(Trajectory(:,:,ll)-StateMean)+MeanPC2;
    plot(TrajectoryPC1,TrajectoryPC2,'LineWidth',1);
end
hold off

%% Estimate preferred orientation

C = 4;
epsilon = 0.2;
mu = 0.3*pi;

NInputSample = 50;
% dSample = 2*pi/NInputSample;
% SampleInput = 0:dSample:2*pi;
NInputSample = NInputSample + 1;
I0 = ExternalInput(NInputSample,Ns,C,epsilon,mu);
%% Euler integration

dt = 1e-3; % s
tmax = 5;
StimTime = .5;
step = round(tmax/dt);

NInputSample = 50;
dSample = 2*pi/NInputSample;
NInputSample = NInputSample + 1;
SampleInput = 0:dSample:2*pi;

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
    MemoryInput = M*MM_old + IEc + MForward*MS_old;
    mMemory_new = qm(MemoryInput);

    SS_old = SS_new;
    SM_old = SM_new;
    MS_old = MS_new;
    MM_old = MM_new;
    mSensory_old = mSensory_new;
    mMemory_old = mMemory_new;
end

PF = ComputePF(mMemory_new,SampleInput);
%%
DecodedTrajectory = zeros(length(PlotTime),5);
for kk = 1:5
    for jj = 1:length(PlotTime)
        DecodedTrajectory(jj,kk) = PVDecoder(PF,MemoryFiringRate(:,jj,kk));
    end
end
DecodedTrajectory = DecodedTrajectory/pi/2*180;
f2 = figure;
figure(f2)
plot(PlotTime',DecodedTrajectory,'LineWidth',1.5);
box off
set(gca,'FontSize',10,'LooseInset',[0,0,0,0],'LineWidth',0.8,'TickDir','out','TickLength',[0.02,0.01]);
xlabel('Time (s)');
ylabel('\theta (deg)');
xlim([0 20]);
ylim([110,130]);
set(gcf,'Units','Centimeters','Position',[2,2,5,5]);