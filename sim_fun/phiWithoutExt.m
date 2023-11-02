function out = phiWithoutExt(W,s,IEc,qs,qm)
N = length(IEc);
b = zeros(N,1);
b = [b;b;IEc;IEc];
Input = W*s+b;
Upper2N = qs(Input(1:2*N,:));
Lower2N = qm(Input((2*N+1):end,:));
out = [Upper2N;Lower2N];
end