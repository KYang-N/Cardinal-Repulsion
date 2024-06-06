%% Search Jf and Jb to investigate the bias and SD profile for the two-layer model

clear
SaveFlag = 0;
LoadFlag = 1;
%% Define neurons

if LoadFlag == 0
    seed = 439;
    rng(seed);
    Tau_syn = .01;
    Ns = 300;
    Nm = 300;
    ParallelNum = 16;
    NE = 2; th = 0.1; sig = 6; maxf = 100;
    qs = @(x) maxf*(x-th).^NE./(sig^NE+(x-th).^NE).*(x>th);
%%  Generate sensory circuit connectivity

    JI = 0.35;
    JE = 0.6;
    lambda = 0.36*pi;
    lambdaI = 1.1*pi;
    dthetas = 2*pi/Ns;
    thetas = 0:dthetas:2*pi-dthetas;
    alpha = 0.04;
    Conn = zeros(Ns,Ns);

    for i = 1:Ns
        thetaCurrent = thetas(i);
        JEModulated = JE*ConnModulation(thetaCurrent,-alpha);
        J = (-JI + JEModulated*exp(-((thetas-pi)/lambda).^2))/(2*pi);
        Conn(i,:) = dthetas*circshift(J,[0 round(-pi/dthetas-1+i)]);
    end
%% Define memory neurons

    dthetam = 2*pi/Nm;
    thetam = 0:dthetam:2*pi-dthetam;
    % Input-output transfer function
    NEM = 1.5; thM = 0.1; sigM = 6.6; maxf = 100;
    qm = @(x) maxf*(x-thM).^NEM./(sigM^NEM+(x-thM).^NEM).*(x>thM);
%% Generate memory circuit connectivity

    JEM = 1; JIM = 0.17;
    lambdaM = 0.2*pi;

    ME = zeros(Nm,Nm);
    MI = ME;
    fE = JEM*exp(-((thetam-pi)/lambdaM).^2)/(2*pi);
    fI = JIM*exp(-((thetam-pi)/(pi*0.6)).^2)/(2*pi);
    for i = 1:Nm
        ME(i,:) = dthetam*circshift(fE,[0 -round(pi/dthetam)-1+i]);
        MI(i,:) = dthetam*circshift(fI,[0 -round(pi/dthetam)-1+i]);
    end

    M = ME-MI;

    IEc = 0.6*ones(Nm,1);
%% Define connectivity and background input

    JfRange = linspace(0.1,0.4,21);
    JbRange = linspace(0.1,0.4,21);

    C = 4;
    epsilon = 0.2;
    mu = 0.3*pi;

    dt = 1e-3; % s
    tmax = 1.5;
    RepTime = 3e3;
    InputOrientation = [0 pi/4 pi/2];
    NumInput = length(InputOrientation);
    step = round(tmax/dt);
    SDIndex = zeros(length(JfRange),length(JbRange));
    MaxBias = SDIndex;
    DecodedOrientation = zeros(RepTime,NumInput);
    tic
    p = parpool(ParallelNum);
    for j = 1:length(JfRange)
        for i = 1:length(JbRange)
            % Generate feedforward connectivity
            JForward = JfRange(j);
            nu = 0.17*pi;
            MForward = zeros(Nm,Ns);
            fForward = 1/(2*pi)*JForward*exp(-((pi-thetas)/nu).^2);
            for ii= 1:Nm
                MForward(ii,:) = dthetas*circshift(fForward,[0 ...
                    -round(pi/dthetas)-1+ii]);
            end

            % Generate feedback connectivity
            JBackward = JbRange(i);
            nufdbk = 0.17*pi;
            dthetam = 2*pi/Nm;
            thetam = 0:dthetam:2*pi-dthetam;
            MBackward = zeros(Ns,Nm);
            fBackward = 1/(2*pi)*JBackward*exp(-((pi-thetam)/nufdbk).^2);
            for jj= 1:Ns
                MBackward(jj,:) = dthetam*circshift(fBackward,[0 ...
                    -round(pi/dthetam)-1+jj]);
            end

            % Compute preferred orietations
            NInput = 50;
            dSample = 2*pi/NInput;
            SampleInput = 0:dSample:2*pi;
            NInput = NInput + 1;

            PFTime = 5; % Saturation
            StimTime = PFTime;
            I0 = ExternalInput(NInput,Ns,C,epsilon,mu);

            % Synaptic variables
            SS_old = zeros(Ns,NInput);
            SM_old = zeros(Nm,NInput);
            MS_old = zeros(Ns,NInput);
            MM_old = zeros(Nm,NInput);
            % Firing rates
            mSensory_old = zeros(Ns,NInput);
            mMemory_old = zeros(Nm,NInput);
            for ii = 1:round(PFTime/dt)
                SS_new = SS_old+ 1/Tau_syn*dt*(-SS_old+mSensory_old);
                SM_new = SM_old + 1/Tau_syn*dt*(-SM_old+mMemory_old);
                MS_new = MS_old + 1/Tau_syn*dt*(-MS_old+mSensory_old);
                MM_new = MM_old + 1/Tau_syn*dt*(-MM_old+mMemory_old);
                SensoryInputI = Conn*SS_old + MBackward*SM_old+I0'*(ii<(StimTime/dt));

                mSensory_new = qs(SensoryInputI);
                MemoryInput = M*MM_old + IEc + MForward*MS_old;
                mMemory_new1 = qm(MemoryInput);
                SS_old = SS_new;
                SM_old = SM_new;
                MS_old = MS_new;
                MM_old = MM_new;
                mSensory_old = mSensory_new;
                mMemory_old = mMemory_new1;
            end
            PF = ComputePF(mMemory_new1,SampleInput);
            StimTime = .5;

            NInput = length(InputOrientation);
            parfor jj = 1:RepTime
                % Synaptic variables
                SS_old = zeros(Ns,NInput);
                SM_old = zeros(Nm,NInput);
                MS_old = zeros(Ns,NInput);
                MM_old = zeros(Nm,NInput);
                % Firing rates
                mSensory_old = zeros(Ns,NInput);
                mMemory_old = zeros(Nm,NInput);
                I0 = zeros(NInput,Ns);
                for kk = 1:NInput
                    I_ext = C*(1-2*epsilon+2*epsilon*exp(-((thetas-pi)/mu).^2));
                    I0(kk,:) = circshift(I_ext,round(((-pi+InputOrientation(kk))/dthetas)));
                end
                % Euler-Maruyama integration
                for ii = 1:step
                    SS_new = SS_old+ 1/Tau_syn*dt*(-SS_old+mSensory_old)+...
                        1/Tau_syn*sqrt(dt)*sqrt(mSensory_old*dt).*randn(Ns,NumInput); 
                    SM_new = SM_old + 1/Tau_syn*dt*(-SM_old+mMemory_old)+...
                        1/Tau_syn*sqrt(dt)*sqrt(mMemory_old*dt).*randn(Nm,NumInput);
                    MS_new = MS_old + 1/Tau_syn*dt*(-MS_old+mSensory_old)+...
                        1/Tau_syn*sqrt(dt)*sqrt(mSensory_old*dt).*randn(Ns,NumInput);
                    MM_new = MM_old + 1/Tau_syn*dt*(-MM_old+mMemory_old)+...
                        1/Tau_syn*sqrt(dt)*sqrt(mMemory_old*dt).*randn(Nm,NumInput);
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
                for kk = 1:NumInput
                    DecodedOrientation(jj,kk) = PVDecoder(PF,mMemory_new(:,kk));
                end
            end
            Bias = DecodedOrientation-InputOrientation;
            Bias(Bias<-pi) = Bias(Bias<-pi) + 2*pi;
            Bias(Bias>pi) = Bias(Bias>pi) - 2*pi;
            Bias = Bias/pi*180/2;
            SD = std(Bias,0,1);
            MaxBias(j,i) = mean(Bias(:,2));
            SDIndex(j,i) = (SD(end) - SD(1))/(SD(end)+SD(1));
        end

        disp([num2str(j),' Jf value(s) finished.'])
        toc
    end
    delete(p)
    if SaveFlag 
        save('TwoLayerBiasSDSearchJfJb.mat');
    end
end
%% Visualize results

if LoadFlag
    DataDir = '';
    load([DataDir,'TwoLayerBiasSDSearchJfJb.mat']);
end

f1 = figure;
figure(f1)
imagesc(JbRange,JfRange,MaxBias);
h = colorbar;
set(h,'LineWidth',0.8,'Ticks',[0 4 8]);
h.TickDirection = 'out';
xlabel('J_b','FontSize',10);
ylabel('J_f','FontSize',10);
set(gca,'Ydir','normal','FontSize',10,'LineWidth',0.8,'TickDir','out', ...
    'TickLength',[0.02,0.01],'LooseInset',[0 0 0 0]);
yticks(0.1:0.15:0.4);
xticks(0.1:0.15:0.4);
clim([0,8]);
box off
axis square
set(gcf,'Unit','Centimeters','Position',[2,2,5,5]);

f2 = figure;
figure(f2)
imagesc(JbRange,JfRange,SDIndex);
h = colorbar;
h.TickDirection = 'out';
h.Ticks = [-0.5,0,0.2];
set(h,'LineWidth',0.8);
clim([-0.5 0.2]);
set(gca,'Ydir','normal','FontSize',10,'LineWidth',0.8,'TickDir', ...
    'out','TickLength',[0.02,0.01],'LooseInset',[0 0 0 0]);
xlabel('J_b','FontSize',10);
ylabel('J_f','FontSize',10);
yticks(0.1:0.15:0.4);
xticks(0.1:0.15:0.4);
set(gcf,'Unit','Centimeters','Position',[2,2,5,5]);
box off
axis square