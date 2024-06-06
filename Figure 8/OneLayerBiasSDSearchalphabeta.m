%% Search alpha and beta to investigate the bias and SD patterns for the one-layer model

clear
SaveFlag = 0;
LoadFlag = 1;
%% Define neurons

if LoadFlag == 0
    seed = 439;
    rng(seed);
    Nm = 300;
    ParallelNum = 4; % Number of parallel workers
    dthetam = 2*pi/Nm;
    thetam = 0:dthetam:2*pi-dthetam;
    Tau_syn = .01;
    % Input-output transfer function
    NEM = 1.5; thM = 0.1; sigM = 6.6; maxf = 100;
    qm = @(x) maxf*(x-thM).^NEM./(sigM^NEM+(x-thM).^NEM).*(x>thM);
%% Define connectivity and background input

    alphaRange = linspace(1e-4,1e-3,21);
    betaRange = linspace(5e-4,5e-3,21);

    JE = 1; JI = .17;
    lambdaM = 0.2*pi;

    dt = 1e-3; % s
    tmax = 1.5;
    RepTime = 3e3;
    InputOrientation = [0,pi/4,pi/2];
    NumInput = length(InputOrientation);
    IEc = 0.6*ones(Nm,1);
    step = round(tmax/dt);
    SDIndex = zeros(length(alphaRange),length(betaRange));
    MaxBias = SDIndex;
    DecodedOrientation = zeros(RepTime,NumInput);

    ME = zeros(Nm,Nm);
    MI = zeros(Nm,Nm);
 
    tic
    p = parpool(ParallelNum);
    for j = 1:length(alphaRange)
        for i = 1:length(betaRange)
            
            % Generate sensory circuit connectivity
            alpha = alphaRange(j);
            beta = betaRange(i);

            for ii = 1:Nm
                thetaCurrent = thetam(ii);
                JEModulated = JE*ConnModulation(thetaCurrent,alpha);
                JIModulated = JI*ConnModulation(thetaCurrent,beta);
                fE = JEModulated*exp(-((thetam-pi)/lambdaM).^2)/(2*pi);
                fI = JIModulated*exp(-((thetam-pi)/(pi*0.6)).^2)/(2*pi);
                ME(ii,:) = dthetam*circshift(fE,[0 round(-pi/dthetam-1+ii)]);
                MI(ii,:) = dthetam*circshift(fI,[0 -round(pi/(dthetam))-1+ii]);
            end
            M = ME - MI;

            % Compute preferred orietations
            NInput = 50;
            dSample = 2*pi/NInput;
            SampleInput = 0:dSample:2*pi;
            NInput = NInput + 1;

            PFTime = 5; % To reach steady state
            StimTime = PFTime;
            mMemory_old = zeros(Nm,NInput);
            MM_old = zeros(Nm,NInput);
            I0 = zeros(NInput,Nm);
            for jj = 1:NInput
                Iext = (cos(thetam-pi)+1)/2;
                I0(jj,:)= circshift(Iext,round(((-pi+SampleInput(jj))/dthetam)));
            end
            for ii = 1:PFTime
                MM_new = MM_old + 1/Tau_syn*dt*(-MM_old+mMemory_old);
                MemoryInput = M*MM_old + IEc + I0'*(ii<(StimTime/dt));
                mMemory_new1 = qm(MemoryInput);
                mMemory_old = mMemory_new1;
                MM_old = MM_new;
            end
            PF = ComputePF(mMemory_new1,SampleInput);
            StimTime = .5;

            parfor jj = 1:RepTime
                mMemory_old = zeros(Nm,NumInput);
                MM_old = zeros(Nm,NumInput);
                mMemory_new = zeros(Nm,NumInput);
                I0 = zeros(NumInput,Nm);
                % Input
                for kk = 1:NumInput
                    I_ext = (cos(thetam-pi)+1)/2;
                    I0(kk,:) = circshift(I_ext,round(((-pi+InputOrientation(kk))/dthetam)));
                end
                % Euler-Maruyama integration
                for ii = 1:step
                    MM_new = MM_old + 1/Tau_syn*dt*(-MM_old+mMemory_old)+  ...
                        1/Tau_syn*sqrt(dt)*sqrt(mMemory_old*dt).*randn(Nm,NumInput);
                    MemoryInput = M*MM_old + IEc + I0'*(ii<(StimTime/dt));
                    mMemory_new = qm(MemoryInput);
                    mMemory_old = mMemory_new;
                    MM_old = MM_new;
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

        disp([num2str(j),' alpha value(s) finished.'])
        toc
    end
    delete(p)
    if SaveFlag 
        save('OneLayerBiasSDSearchalphabeta.mat');
    end
end
%% Visualize results

if LoadFlag
    clear %#ok<*UNRCH> 
    load OneLayerBiasSDSearchalphabeta.mat;
end

f1 = figure;
figure(f1)
imagesc(betaRange,alphaRange,MaxBias);
h = colorbar;
set(h,'LineWidth',0.8);
h.TickDirection = 'out';
set(gca,'Ydir','normal','FontSize',10,'LineWidth',0.8,'TickDir','out', ...
    'TickLength',[0.025,0.01],'LooseInset',[0 0 0 0]);
xlabel('\beta','FontSize',10);
ylabel('\alpha','FontSize',10);
yticks(1e-4:4.5e-4:1e-3);
xticks(5e-4:2.25e-3:5e-3);
clim([-4,4]);
h.Ticks = [-4,0,4];
box off
axis square
set(gcf,'Unit','Centimeters','Position',[2,2,5,5]);

f2 = figure;
figure(f2)
imagesc(betaRange,alphaRange,SDIndex);
h = colorbar;
h.TickDirection = 'out';
h.Ticks = [-0.15 0 0.1];
set(h,'LineWidth',0.8);
clim([-0.15 0.1]);
set(gca,'Ydir','normal','FontSize',10,'LineWidth',0.8,'TickDir', ...
    'out','TickLength',[0.025,0.01],'LooseInset',[0 0 0 0]);
xlabel('\beta','FontSize',10);
ylabel('\alpha','FontSize',10);
yticks(1e-4:4.5e-4:1e-3);
xticks(5e-4:2.25e-3:5e-3);
box off
axis square
set(gcf,'Unit','Centimeters','Position',[2,2,5,5]);