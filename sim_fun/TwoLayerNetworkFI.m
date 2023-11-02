function [FI,SensoryNet,MemoryNet,DynParams] = TwoLayerNetworkFI(SensoryNet,MemoryNet,DynParams)

if ~isfield(DynParams,'NoiseTime')
    DynParams.NoiseTime = 0;
end

if ~isfield(DynParams,'Parallel')
    DynParams.Parallel = 0;
end

if ~isfield(DynParams,'AddNoise')
    DynParams.AddNoise = 1;
end

if ~isfield(DynParams,'RepTime')
    DynParams.RepTime = 1e3;
end

if ~isfield(SensoryNet,'tau')
    SensoryNet.tau = 10e-3;
end

if ~isfield(SensoryNet,'labmda')
    SensoryNet.lambda = 0.36*pi;
end

if ~isfield(SensoryNet,'lambdaI')
    SensoryNet.lambdaI = 1.1*pi;
end

if ~isfield(DynParams,'NInputSample')
    DynParams.NInputSample = round(SensoryNet.N/6) + 1;
end

if ~isfield(DynParams,'dSample')
    DynParams.dSample = 2*pi/(DynParams.NInputSample-1);
end

if ~isfield(DynParams,'C')
    DynParams.C = 4;
end

if ~isfield(DynParams,'mu')
    DynParams.mu = 0.3*pi;
end

if ~isfield(DynParams,'epsilon')
    DynParams.epsilon = 0.2;
end

if ~isfield(SensoryNet,'q')
    NE = 2; th = 0.1; sig = 6; maxf = 100;
    SensoryNet.NE = NE; SensoryNet.th = th; SensoryNet.sig = sig; SensoryNet.maxf = maxf;
    SensoryNet.q = @(x) maxf*(x-th).^NE./(sig^NE+(x-th).^NE).*(x>th);
end

if ~isfield(MemoryNet,'q')
    NEM = 1.5; thM = 0.1; sigM = 6.6; maxf = 100;
    MemoryNet.thM= thM; MemoryNet.NEM = NEM; MemoryNet.sigM = sigM;
    MemoryNet.q = @(x) maxf*(x-thM).^NEM./(sigM.^NEM+(x-thM).^NEM).*(x>thM);
end

if ~isfield(SensoryNet,'JI')
    SensoryNet.JI = .35;
end

if ~isfield(SensoryNet,'JE')
    SensoryNet.JE = .60;
end

if ~isfield(SensoryNet,'Conn')
    SensoryNet = SensoryNetRecurConn(SensoryNet);
end

if ~isfield(MemoryNet,'JE')
    MemoryNet.JE = 1;
end

if ~isfield(MemoryNet,'JI')
    MemoryNet.JI = 0.17;
end

if ~isfield(MemoryNet,'lambdaM')
    MemoryNet.lambdaM = 0.2*pi;
end

if ~isfield(MemoryNet,'Conn')
    MemoryNet = MemoryNetRecurConn(MemoryNet);
end

if ~isfield(MemoryNet,'IEc')
    MemoryNet.IEc = 0.6*ones(MemoryNet.N,1);
end

if ~isfield(DynParams,'FullModel_tmax')
    DynParams.FullModel_tmax = max(DynParams.FullDecodeTime);
end

Tau_syn = SensoryNet.tau;
NInputSample = DynParams.NInputSample;
SampleInput = 0:DynParams.dSample:2*pi;
Ns = SensoryNet.N;
Nm = MemoryNet.N;
% dSample = DynParams.dSample;
C = DynParams.C;
epsilon = DynParams.epsilon;
mu = DynParams.mu;
tmax = DynParams.FullModel_tmax;
dt = DynParams.dt;
StimTime = DynParams.StimTime;
DecodeTime = DynParams.FullDecodeTime;

TuningCurves1 = zeros(MemoryNet.N,NInputSample,length(DecodeTime));
TuningCurves2 = zeros(MemoryNet.N,NInputSample,length(DecodeTime));
TuningCurves3 = zeros(MemoryNet.N,NInputSample,length(DecodeTime));

I0 = ExternalInput(NInputSample,SensoryNet.N,C,epsilon,mu);

step = round(tmax/dt);

%% SDE Dynamics
% Synaptic variables
Start = datetime("now");
if DynParams.Parallel == 0
    for ll = 1:DynParams.RepTime
        DecodeCounter = 1;
        SS_old = zeros(SensoryNet.N,NInputSample);
        SM_old = zeros(MemoryNet.N,NInputSample);
        MS_old = zeros(SensoryNet.N,NInputSample);
        MM_old = zeros(MemoryNet.N,NInputSample);
        % Firing rates
        mSensory_old = zeros(SensoryNet.N,NInputSample);
        mMemory_old = zeros(MemoryNet.N,NInputSample);

        for ii = 1:step
            NoiseFlag = (ii>(DynParams.NoiseTime/dt))*DynParams.AddNoise;
            SS_new = SS_old + 1/Tau_syn*dt*(-SS_old+mSensory_old)+...
                1/Tau_syn*sqrt(dt)*sqrt(mSensory_old*dt).*randn(Ns,NInputSample)*NoiseFlag;
            SM_new = SM_old + 1/Tau_syn*dt*(-SM_old+mMemory_old)+...
                1/Tau_syn*sqrt(dt)*sqrt(mMemory_old*dt).*randn(Nm,NInputSample)*NoiseFlag;
            MS_new = MS_old + 1/Tau_syn*dt*(-MS_old+mSensory_old)+...
                1/Tau_syn*sqrt(dt)*sqrt(mSensory_old*dt).*randn(Ns,NInputSample)*NoiseFlag;
            MM_new = MM_old + 1/Tau_syn*dt*(-MM_old+mMemory_old)+...
                1/Tau_syn*sqrt(dt)*sqrt(mMemory_old*dt).*randn(Nm,NInputSample)*NoiseFlag;
            SensoryInputI = SensoryNet.Conn*SS_old + MemoryNet.MBackward*SM_old ...
                + I0'*(ii<(StimTime/dt));
            mSensory_new = SensoryNet.q(SensoryInputI);
            MemoryInput = MemoryNet.Conn*MM_old + MemoryNet.IEc + SensoryNet.MForward*MS_old;
            mMemory_new = MemoryNet.q(MemoryInput);

            if ii == round(DecodeTime(DecodeCounter)/dt)
                switch DecodeCounter
                    case 1
                        TuningCurves1(:,:,ll) = mMemory_new;
                    case 2
                        TuningCurves2(:,:,ll) = mMemory_new;
                    case 3
                        TuningCurves3(:,:,ll) = mMemory_new;
                end
                DecodeCounter = DecodeCounter + 1;
            end
            SS_old = SS_new;
            SM_old = SM_new;
            MS_old = MS_new;
            MM_old = MM_new;
            mSensory_old = mSensory_new;
            mMemory_old = mMemory_new;
        end
    end

elseif DynParams.Parallel ~= 0
    parfor ll = 1:DynParams.RepTime
        DecodeCounter = 1;
        SS_old = zeros(SensoryNet.N,NInputSample);
        SM_old = zeros(MemoryNet.N,NInputSample);
        MS_old = zeros(SensoryNet.N,NInputSample);
        MM_old = zeros(MemoryNet.N,NInputSample);
        % Firing rates
        mSensory_old = zeros(SensoryNet.N,NInputSample);
        mMemory_old = zeros(MemoryNet.N,NInputSample);

        for ii = 1:step
            NoiseFlag = (ii>(DynParams.NoiseTime/dt))*DynParams.AddNoise;
            SS_new = SS_old + 1/Tau_syn*dt*(-SS_old+mSensory_old)+...
                1/Tau_syn*sqrt(dt)*sqrt(mSensory_old*dt).*randn(Ns,NInputSample)*NoiseFlag;
            SM_new = SM_old + 1/Tau_syn*dt*(-SM_old+mMemory_old)+...
                1/Tau_syn*sqrt(dt)*sqrt(mMemory_old*dt).*randn(Nm,NInputSample)*NoiseFlag;
            MS_new = MS_old + 1/Tau_syn*dt*(-MS_old+mSensory_old)+...
                1/Tau_syn*sqrt(dt)*sqrt(mSensory_old*dt).*randn(Ns,NInputSample)*NoiseFlag;
            MM_new = MM_old + 1/Tau_syn*dt*(-MM_old+mMemory_old)+...
                1/Tau_syn*sqrt(dt)*sqrt(mMemory_old*dt).*randn(Nm,NInputSample)*NoiseFlag;
            SensoryInputI = SensoryNet.Conn*SS_old + MemoryNet.MBackward*SM_old ...
                + I0'*(ii<(StimTime/dt));
            mSensory_new = SensoryNet.q(SensoryInputI);
            MemoryInput = MemoryNet.Conn*MM_old + MemoryNet.IEc + SensoryNet.MForward*MS_old;
            mMemory_new = MemoryNet.q(MemoryInput);

            if ii == round(DecodeTime(DecodeCounter)/dt)
                switch DecodeCounter
                    case 1
                        TuningCurves1(:,:,ll) = mMemory_new;
                    case 2
                        TuningCurves2(:,:,ll) = mMemory_new;
                    case 3
                        TuningCurves3(:,:,ll) = mMemory_new;
                end
                DecodeCounter = DecodeCounter + 1;
            end
            SS_old = SS_new;
            SM_old = SM_new;
            MS_old = MS_new;
            MM_old = MM_new;
            mSensory_old = mSensory_new;
            mMemory_old = mMemory_new;
        end
    end
end
FI = zeros(3,NInputSample);
for uu = 1:3
    switch uu
        case 1
            TCMean = mean(TuningCurves1,3);
            TCVar = var(TuningCurves1,0,3);
        case 2
            TCMean = mean(TuningCurves2,3);
            TCVar = var(TuningCurves2,0,3);
        case 3
            TCMean = mean(TuningCurves3,3);
            TCVar = var(TuningCurves3,0,3);
    end
    FI(uu,:) = ComputeFI(SampleInput,TCMean,TCVar);
end
Lap = datetime("now");
disp(['Time elapsed: ',datestr(Lap-Start,'HH:MM:SS')]) %#ok<*DATST>
end