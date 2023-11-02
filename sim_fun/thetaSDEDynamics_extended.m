function DecodedOrientation = thetaSDEDynamics_extended(av,bv,aD,bD,alambda,blambda,NFourier,DynParams,ThetaLoc)
% lambda is the largest left eigenvalues of the Jacobian matrix K
DecodeTime = DynParams.SimplifiedDecodeTime;
DecodedOrientation = zeros(DynParams.RepTime,DynParams.NInputSample, ...
    length(DecodeTime));
step = round(DynParams.Manifold_tmax/DynParams.dt);
Start = datetime("now");
dt = DynParams.dt;
for ll = 1:DynParams.RepTime
    DecodeCounter = 1;
    theta_old = ThetaLoc;
    delta_theta = 0;
    for ii = 1:step
        theta_new = theta_old+dt*InvFourier(av,bv,theta_old,NFourier)+...
             sqrt(dt)*sqrt(2*InvFourier(aD,bD,theta_old,NFourier)).*...
        randn(1,DynParams.NInputSample) + dt*InvFourier(alambda,blambda,theta_old,NFourier).*(delta_theta);

        if ii == round(DecodeTime(DecodeCounter)/DynParams.dt)
            DecodedOrientation(ll,:,DecodeCounter) = theta_new;
            if DecodeCounter < length(DecodeTime)
                DecodeCounter = DecodeCounter + 1;
            end
        end
        delta_theta = theta_new-theta_old;
        theta_old = theta_new;
    end
end
Lap = datetime("now");
disp(['Time elapsed: ',datestr(Lap-Start,'HH:MM:SS')]) %#ok<*DATST>