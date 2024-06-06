function dsdthetaNorm = ComputeTVNormOneLayer(MemoryNet,DynParams)
Tau_syn = MemoryNet.tau;
NInputSample = DynParams.NInputSample;
dSample = DynParams.dSample;
SampleInput = 0:dSample:2*pi;

tmax = DynParams.Manifold_tmax;
dt = DynParams.dt;
StimTime = DynParams.StimTime;

I0 = ExternalInput_OneLayer(NInputSample,MemoryNet.N);

% Synaptic variables
MM_old = zeros(MemoryNet.N,NInputSample);
% Firing rates
mMemory_old = zeros(MemoryNet.N,NInputSample);
step = round(tmax/dt);

for ii = 1:step
    MM_new = MM_old + 1/Tau_syn*dt*(-MM_old+mMemory_old);
    MemoryInput = MemoryNet.Conn*MM_old + MemoryNet.IEc + I0'*(ii<(StimTime/dt));
    mMemory_new = MemoryNet.q(MemoryInput);
    MM_old = MM_new;
    mMemory_old = mMemory_new;
end

PF = ComputePF(mMemory_old,SampleInput);

% ThetaLoc = zeros(1,NInputSample);
% for ii = 1:NInputSample
%     ThetaLoc(ii) = PVDecoder(PF,mMemory_new(:,ii));
% end

s_manifold = MM_old;

SampleInput = SampleInput/2*180/pi; % Convert to degrees
sRight = [s_manifold(:,2:end),s_manifold(:,2)]; % note that s(:,1) == s(:,end)
sLeft = [s_manifold(:,end-1),s_manifold(:,1:end-1)];
ds = sRight-sLeft;
ThetaRight = [SampleInput(2:end),SampleInput(2)];
ThetaLeft = [SampleInput(end-1),SampleInput(1:end-1)];
dtheta = ThetaRight-ThetaLeft;
dtheta(dtheta<0) = dtheta(dtheta<0) + 180;
deri = ds./dtheta;

dsdthetaNorm = zeros(1,NInputSample);
for i = 1:NInputSample
    dsdthetaNorm(i) = norm(deri(:,i));  % s
end


