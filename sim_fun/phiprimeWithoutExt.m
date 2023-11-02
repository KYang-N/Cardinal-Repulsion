function out = phiprimeWithoutExt(W,s,SensoryNet,MemoryNet)
th = MemoryNet.thM;
maxf = SensoryNet.maxf;
NE = SensoryNet.NE;
NEM = MemoryNet.NEM;
sig = SensoryNet.sig;
sigM = MemoryNet.sigM;
qsprime = @(x) NE*maxf*(x-th).^(NE-1)*sig^NE./((sig^NE+(x-th).^NE).^2).*(x>th);
qmprime = @(x) NEM*maxf*((x-th).*(x>th)).^(NEM-1)*sigM^NEM./((sigM^NEM+((x-th).*(x>th)).^NEM).^2);
N = length(MemoryNet.IEc);
b = zeros(N,1);
b = [b;b;MemoryNet.IEc;MemoryNet.IEc];
Input = W*s+b;
Upper2N = qsprime(Input(1:2*N));
Lower2N = qmprime(Input((2*N+1):end));
out = [Upper2N;Lower2N];
end