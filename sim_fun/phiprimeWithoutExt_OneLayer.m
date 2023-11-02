function out = phiprimeWithoutExt_OneLayer(W,s,MemoryNet)
th = MemoryNet.thM;
maxf = MemoryNet.maxf;
NEM = MemoryNet.NEM;
sigM = MemoryNet.sigM;
qmprime = @(x) NEM*maxf*sqrt((x-th).*(x>th))*sigM^NEM./(sigM^NEM+((x-th).*(x>th)).^NEM).^2.*(x>th);
Input = W*s+MemoryNet.IEc;
out = qmprime(Input);
end