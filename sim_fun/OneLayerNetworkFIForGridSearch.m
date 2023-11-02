function [FI,MemoryNet,DynParams] = OneLayerNetworkFIForGridSearch(MemoryNet,DynParams)

if ~isfield(DynParams,'NoiseTime')
    DynParams.NoiseTime = 0;
end

if ~isfield(DynParams,'AddNoise')
    DynParams.AddNoise = 1;
end

if ~isfield(DynParams,'RepTime')
    DynParams.RepTime = 1;
end

if ~isfield(MemoryNet,'tau')
    MemoryNet.tau = 1e-2;
end

if ~isfield(DynParams,'NInputSample')
    DynParams.NInputSample = round(MemoryNet.N/6) + 1;
end

if ~isfield(DynParams,'dSample')
    DynParams.dSample = 2*pi/(DynParams.NInputSample-1);
end

if ~isfield(MemoryNet,'q')
    NEM = 1.5; thM = 0.1; sigM = 6.6; maxf = 100;
    MemoryNet.thM= thM; MemoryNet.NEM = NEM; MemoryNet.sigM = sigM;
    MemoryNet.q = @(x) maxf*(x-thM).^NEM./(sigM.^NEM+(x-thM).^NEM).*(x>thM);
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
    MemoryNet = OneLayerRecurConn(MemoryNet);
end

if ~isfield(MemoryNet,'IEc')
    MemoryNet.IEc = 0.6*ones(MemoryNet.N,1);
end


Tau_syn = MemoryNet.tau;
NInputSample = DynParams.NInputSample;
SampleInput = 0:DynParams.dSample:2*pi;
Nm = MemoryNet.N;
dthetam = 2*pi/Nm;
thetam = 0:dthetam:2*pi-dthetam;
% dSample = DynParams.dSample;
tmax = DynParams.DecodeTime;
dt = DynParams.dt;
StimTime = DynParams.StimTime;
DecodeTime = DynParams.DecodeTime;

TuningCurves = zeros(MemoryNet.N,NInputSample,DynParams.RepTime);


I0 = zeros(NInputSample,Nm);
for jj = 1:NInputSample
    Iext = (cos(thetam-pi)+1)/2;
    I0(jj,:) = circshift(Iext,round(((-pi+SampleInput(jj))/dthetam)));
end
step = round(tmax/dt);
%% SDE Dynamics
% Synaptic variables


if DynParams.Parallel == 0
    Start = datetime("now");
    for ll = 1:DynParams.RepTime
        MM_old = zeros(Nm,NInputSample);
        % Firing rates
        mMemory_old = zeros(Nm,NInputSample);

        for ii = 1:step
            NoiseFlag = (ii>(DynParams.NoiseTime/dt))*DynParams.AddNoise;
            MM_new = MM_old + 1/Tau_syn*dt*(-MM_old+mMemory_old)+...
                1/Tau_syn*sqrt(dt)*sqrt(mMemory_old*dt).*randn(Nm,NInputSample)*NoiseFlag;
            MemoryInput = MemoryNet.Conn*MM_old + MemoryNet.IEc + I0'*(ii<(StimTime/dt));
            mMemory_new = MemoryNet.q(MemoryInput);

            if ii == round(DecodeTime/dt)
                TuningCurves(:,:,ll) = mMemory_new;

            end
            MM_old = MM_new;
            mMemory_old = mMemory_new;
        end
    end

else

    Start = datetime("now");
    parfor ll = 1:DynParams.RepTime
        MM_old = zeros(Nm,NInputSample);
        % Firing rates
        mMemory_old = zeros(Nm,NInputSample);

        for ii = 1:step
            NoiseFlag = (ii>(DynParams.NoiseTime/dt))*DynParams.AddNoise;
            MM_new = MM_old + 1/Tau_syn*dt*(-MM_old+mMemory_old)+...
                1/Tau_syn*sqrt(dt)*sqrt(mMemory_old*dt).*randn(Nm,NInputSample)*NoiseFlag;
            MemoryInput = MemoryNet.Conn*MM_old + MemoryNet.IEc + I0'*(ii<(StimTime/dt));
            mMemory_new = MemoryNet.q(MemoryInput);

            if ii == round(DecodeTime/dt)

                TuningCurves(:,:,ll) = mMemory_new;

            end
            MM_old = MM_new;
            mMemory_old = mMemory_new;
        end
    end
end

TCMean = mean(TuningCurves,3);
TCVar = var(TuningCurves,0,3);

FI = ComputeFI(SampleInput,TCMean,TCVar);

Lap = datetime("now");
disp(['Time elapsed: ',datestr(Lap-Start,'HH:MM:SS')]) %#ok<*DATST>
end