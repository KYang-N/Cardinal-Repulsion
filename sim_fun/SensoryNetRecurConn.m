function SensoryNet = SensoryNetRecurConn(SensoryNet)
dthetas = 2*pi/SensoryNet.N;
thetas = 0:dthetas:2*pi-dthetas;
ConnE = zeros(SensoryNet.N,SensoryNet.N);
ConnI = zeros(SensoryNet.N,SensoryNet.N);

%  Modify gains
if strcmp(SensoryNet.Mode,'EOnly') || strcmp(SensoryNet.Mode,'E')
    for i = 1:SensoryNet.N
        thetaCurrent = thetas(i);
        JEModulated = SensoryNet.JE*ConnModulation(thetaCurrent,-SensoryNet.alpha);
        JIModulated = SensoryNet.JI;
        J = (-JIModulated + JEModulated*exp(-((thetas-pi)/SensoryNet.lambda).^2))/(2*pi);
        SensoryNet.Conn(i,:) = dthetas*circshift(J,[0 round(-pi/dthetas-1+i)]);
    end
elseif strcmp(SensoryNet.Mode,'Both') || strcmp(SensoryNet.Mode,'both')
    for i = 1:SensoryNet.N
        thetaCurrent = thetas(i);
        JEModulated = SensoryNet.JE*ConnModulation(thetaCurrent,SensoryNet.alpha);
        JIModulated = SensoryNet.JI*ConnModulation(thetaCurrent,SensoryNet.beta);
        ConnI(i,:) = 1/2/pi*dthetas*circshift(-JIModulated*exp(-((thetas-pi)/SensoryNet.lambdaI).^2), ...
            [0 round(-pi/dthetas-1+i)]);
        ConnE(i,:) = 1/2/pi*dthetas*circshift(JEModulated*exp(-((thetas-pi)/SensoryNet.lambda).^2), ...
            [0,round(-pi/dthetas-1+i)]);
    end
    SensoryNet.Conn = ConnE+ConnI;
end
