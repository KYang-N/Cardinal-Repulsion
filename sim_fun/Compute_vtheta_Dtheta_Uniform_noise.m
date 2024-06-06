function [vtheta,Dtheta,ThetaLoc] = Compute_vtheta_Dtheta_Uniform_noise(SensoryNet,MemoryNet,DynParams)

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

if ~isfield(DynParams,'DecodingFrom')
    DynParams.DecodingFrom = 'Memory';
end

if ~isfield(DynParams,'sigma')
    DynParams.sigma = 0.2;
end

if ~isfield(SensoryNet,'q')
    NE = 2; th = 0.1; sig = 6; maxf = 100;
    SensoryNet.NE = NE; SensoryNet.th = th; SensoryNet.sig = sig; SensoryNet.maxf = maxf;
    SensoryNet.q = @(x) maxf*(x-th).^NE./(sig^NE+(x-th).^NE).*(x>th);
end

if ~isfield(MemoryNet,'q')
    NEM = 1.5; thM = 0.1; sigM = 6.6; maxf = 100;
    MemoryNet.thM= thM; MemoryNet.NEM = NEM; MemoryNet.sigM = sigM;
    MemoryNet.q = @(x) maxf*(x-thM).^NEM./(sigM^NEM+(x-thM).^NEM).*(x>thM);
end

if ~isfield(SensoryNet,'JI')
    SensoryNet.JI = 0.35;
end

if ~isfield(SensoryNet,'JE')
    SensoryNet.JE = 0.6;
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

Tau_syn = SensoryNet.tau;
NInputSample = DynParams.NInputSample;
dSample = DynParams.dSample;
SampleInput = 0:dSample:2*pi;
C = DynParams.C;
epsilon = DynParams.epsilon;
mu = DynParams.mu;
tmax = DynParams.Manifold_tmax;
dt = DynParams.dt;
StimTime = DynParams.StimTime;

I0 = ExternalInput(NInputSample,SensoryNet.N,C,epsilon,mu);

% Synaptic variables
SS_old = zeros(SensoryNet.N,NInputSample);
SM_old = zeros(MemoryNet.N,NInputSample);
MS_old = zeros(SensoryNet.N,NInputSample);
MM_old = zeros(MemoryNet.N,NInputSample);
% Firing rates
mSensory_old = zeros(SensoryNet.N,NInputSample);
mMemory_old = zeros(MemoryNet.N,NInputSample);
step = round(tmax/dt);

for ii = 1:step
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

if strcmp(DynParams.DecodingFrom,'Memory')
    PF = ComputePF(mMemory_old,SampleInput);
else
    PF = ComputePF(mSensory_old,SampleInput);
end

ThetaLoc = zeros(1,NInputSample);
for ii = 1:NInputSample

    if strcmp(DynParams.DecodingFrom,'Memory')
        ThetaLoc(ii) = PVDecoder(PF,mMemory_new(:,ii));
    else
        ThetaLoc(ii) = PVDecoder(PF,mSensory_new(:,ii));
    end
end

s_manifold = [SS_old;MS_old;MM_old;SM_old];
M_zero = zeros(SensoryNet.N,SensoryNet.N);
W = [SensoryNet.Conn,M_zero,M_zero,MemoryNet.MBackward;SensoryNet.Conn,M_zero,M_zero,MemoryNet.MBackward;
    M_zero,SensoryNet.MForward,MemoryNet.Conn,M_zero;M_zero,SensoryNet.MForward,MemoryNet.Conn,M_zero;];
u = zeros(4*MemoryNet.N,NInputSample);
phiOut = phiWithoutExt(W,s_manifold,MemoryNet.IEc,SensoryNet.q,MemoryNet.q);

dsdtheta = ComputeDerivative(s_manifold,ThetaLoc);
lambda = zeros(1,NInputSample);
for i = 1:NInputSample
    phiPrimeoutput = phiprimeWithoutExt(W,s_manifold(:,i),SensoryNet,MemoryNet);
    K = diag(phiPrimeoutput)*W-eye(4*MemoryNet.N);
    [V,E] = eig(K);
    [lambda(i),Ind] = max(real(diag(E)));
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
    D(i) = 1/(2*(Tau_syn)^2)*sum(uSquare*DynParams.sigma^2);  % s
end

[Dtheta,vtheta] = ConvertTotheta_new(s_manifold,SampleInput,D,v);


