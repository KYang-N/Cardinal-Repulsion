function DecodedOrientation = thetaSDEDynamics(av,bv,aD,bD,NFourier,DynParams,ThetaLoc)
DecodeTime = DynParams.SimplifiedDecodeTime;
DecodedOrientation = zeros(DynParams.RepTime,DynParams.NInputSample, ...
    length(DecodeTime));
dt = DynParams.dt;
step = round(max(DecodeTime)/dt);
Start = datetime("now");

for ll = 1:DynParams.RepTime
    DecodeCounter = 1;
    theta_old = ThetaLoc;
    for ii = 1:step
        theta_new = theta_old+dt*InvFourier(av,bv,theta_old,NFourier)+...
             sqrt(dt)*sqrt(2*InvFourier(aD,bD,theta_old,NFourier)).*...
        randn(1,DynParams.NInputSample);

        if ii == round(DecodeTime(DecodeCounter)/DynParams.dt)
            DecodedOrientation(ll,:,DecodeCounter) = theta_new;
            if DecodeCounter < length(DecodeTime)
                DecodeCounter = DecodeCounter + 1;
            end
        end
        theta_old = theta_new;
    end
end
Lap = datetime("now");
disp(['Time elapsed: ',datestr(Lap-Start,'HH:MM:SS')]) %#ok<*DATST>