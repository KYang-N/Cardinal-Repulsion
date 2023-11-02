function [vtheta,Dtheta] = Compute_vtheta_Dtheta_OneLayer(MemoryNet,DynParams)
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

ThetaLoc = zeros(1,NInputSample);
for ii = 1:NInputSample
    ThetaLoc(ii) = PVDecoder(PF,mMemory_new(:,ii));
end

s_manifold = MM_old;

phiOut = phiWithoutExt_OneLayer(MemoryNet.Conn,s_manifold,MemoryNet.IEc,MemoryNet.q);

dsdtheta = ComputeDerivative(s_manifold,ThetaLoc);
W = MemoryNet.Conn;
u = zeros(MemoryNet.N,NInputSample);

for i = 1:NInputSample
    phiPrimeoutput = phiprimeWithoutExt_OneLayer(W,s_manifold(:,i),MemoryNet);
    K = diag(phiPrimeoutput)*W-eye(MemoryNet.N);
    [V,E] = eig(K);
    [~,Ind] = max(real(diag(E)));
    u(:,i) = V(:,Ind);
    if dsdtheta(:,i)'*u(:,i) < 0
        u(:,i) = -u(:,i);
    end
end

v = zeros(1,NInputSample);
D = zeros(1,NInputSample);
for i = 1:NInputSample
    v(i) = 1/(Tau_syn)*u(:,i)'*(-s_manifold(:,i)+phiOut(:,i)); % s
    uSquare = u(:,i).^2;
    D(i) = 1/(2*(Tau_syn)^2)*uSquare'*phiOut(:,i)*dt;  %s
end

[Dtheta,vtheta] = ConvertTotheta(s_manifold,dSample,D,v);


