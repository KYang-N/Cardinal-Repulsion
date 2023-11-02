function [DecodedOrientation,SensoryNet,MemoryNet,DynParams] = FullSDEDynamicsForComparison(SensoryNet,MemoryNet,DynParams)

if ~isfield(DynParams,'Parallel')
    DynParams.Parallel = 0;
end

if ~isfield(DynParams,'AddNoise')
    DynParams.AddNoise = 1;
end

if ~isfield(DynParams,'RepTime')
    DynParams.RepTime = 3e3;
end

if ~isfield(SensoryNet,'tau')
    SensoryNet.tau = 1e-2;
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

DecodedOrientation = zeros(DynParams.RepTime,NInputSample, ...
    length(DecodeTime));

I0 = ExternalInput(NInputSample,SensoryNet.N,C,epsilon,mu);

%% Estimate preferred orientation
% Synaptic variables
SS_old = zeros(SensoryNet.N,NInputSample);
SM_old = zeros(MemoryNet.N,NInputSample);
MS_old = zeros(SensoryNet.N,NInputSample);
MM_old = zeros(MemoryNet.N,NInputSample);
% Firing rates
mSensory_old = zeros(SensoryNet.N,NInputSample);
mMemory_old = zeros(MemoryNet.N,NInputSample);
step = round(tmax/dt);

for ii = 1:round(DynParams.PFTime/dt)
    SS_new = SS_old + 1/Tau_syn*dt*(-SS_old+mSensory_old);
    SM_new = SM_old + 1/Tau_syn*dt*(-SM_old+mMemory_old);
    MS_new = MS_old + 1/Tau_syn*dt*(-MS_old+mSensory_old);
    MM_new = MM_old + 1/Tau_syn*dt*(-MM_old+mMemory_old);
    SensoryInputI = SensoryNet.Conn*SS_old + MemoryNet.MBackward*SM_old + I0'*(ii<(StimTime/dt));
    mSensory_new = SensoryNet.q(SensoryInputI);
    MemoryInput = MemoryNet.Conn*MM_old + MemoryNet.IEc + SensoryNet.MForward*MS_old;
    mMemory_new = MemoryNet.q(MemoryInput);

    SS_old = SS_new;
    SM_old = SM_new;
    MS_old = MS_new;
    MM_old = MM_new;
    mSensory_old = mSensory_new;
    mMemory_old = mMemory_new;
end
PF = ComputePF(mMemory_new,SampleInput); % preferred orientations

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
                for jj = 1:NInputSample
                    DecodedOrientation(ll,jj,DecodeCounter) = PVDecoder(PF,mMemory_new(:,jj));
                end
                if DecodeCounter < length(DecodeTime)
                    DecodeCounter = DecodeCounter + 1;
                end
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
    p = parpool(DynParams.Parallel);
    parfor ll = 1:DynParams.RepTime
        DecodeCounter = 1;
        DecodedOrientationTemp = zeros(NInputSample,length(DecodeTime));
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
                for jj = 1:NInputSample
                    DecodedOrientationTemp(jj,DecodeCounter) = PVDecoder(PF,mMemory_new(:,jj));
                end
                if DecodeCounter < length(DecodeTime)
                    DecodeCounter = DecodeCounter + 1;
                end
            end
            SS_old = SS_new;
            SM_old = SM_new;
            MS_old = MS_new;
            MM_old = MM_new;
            mSensory_old = mSensory_new;
            mMemory_old = mMemory_new;
        end
        DecodedOrientation(ll,:,:) = DecodedOrientationTemp;
    end
    delete(p);
end
Lap = datetime("now");
disp(['Time elapsed: ',datestr(Lap-Start,'HH:MM:SS')]) %#ok<*DATST>
end