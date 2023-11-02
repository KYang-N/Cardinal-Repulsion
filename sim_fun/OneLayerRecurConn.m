function MemoryNet = OneLayerRecurConn(MemoryNet)
dthetam = 2*pi/MemoryNet.N;
thetam = 0:dthetam:2*pi-dthetam;
ME = zeros(MemoryNet.N,MemoryNet.N);
MI = ME;
if strcmp(MemoryNet.Mode,'EOnly')
    for i = 1:MemoryNet.N
        thetaCurrent = thetam(i);
        JEEModulated = MemoryNet.JE*ConnModulation(thetaCurrent,-MemoryNet.alpha);
        JIModulated = MemoryNet.JI*ConnModulation(thetaCurrent,0);
        fE = JEEModulated*exp(-((thetam-pi)/MemoryNet.lambdaM).^2)/(2*pi);
        fI = JIModulated*exp(-((thetam-pi)/(pi*0.6)).^2)/(2*pi);
        ME(i,:) = dthetam*circshift(fE,[0 -round(pi/dthetam)-1+i]);
        MI(i,:) = dthetam*circshift(fI,[0 -round(pi/dthetam)-1+i]);
    end
elseif strcmp(MemoryNet.Mode,'Both') || strcmp(MemoryNet.Mode,'both')
    for i = 1:MemoryNet.N
        thetaCurrent = thetam(i);
        JEModulated = MemoryNet.JE*ConnModulation(thetaCurrent,MemoryNet.alpha);
        JIModulated = MemoryNet.JI*ConnModulation(thetaCurrent,MemoryNet.beta);
        fE = JEModulated*exp(-((thetam-pi)/MemoryNet.lambdaM).^2)/(2*pi);
        fI = JIModulated*exp(-((thetam-pi)/(pi*0.6)).^2)/(2*pi);
        ME(i,:) = dthetam*circshift(fE,[0 -round(pi/dthetam)-1+i]);
        MI(i,:) = dthetam*circshift(fI,[0 -round(pi/dthetam)-1+i]);
    end
end
MemoryNet.Conn = ME-MI;
end
