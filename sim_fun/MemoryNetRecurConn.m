function MemoryNet = MemoryNetRecurConn(MemoryNet)
dthetam = 2*pi/MemoryNet.N;
thetam = 0:dthetam:2*pi-dthetam;
fE = MemoryNet.JE*exp(-((thetam-pi)/(MemoryNet.lambdaM)).^2)/(2*pi);
fI = MemoryNet.JI*exp(-((thetam-pi)/(pi*0.6)).^2)/(2*pi);

ME = zeros(MemoryNet.N,MemoryNet.N);
MI = ME;

for i = 1:MemoryNet.N
    ME(i,:) = dthetam*circshift(fE,[0 -round(pi/(dthetam))-1+i]);
    MI(i,:) = dthetam*circshift(fI,[0 -round(pi/dthetam)-1+i]);
end

MemoryNet.Conn = ME - MI;
