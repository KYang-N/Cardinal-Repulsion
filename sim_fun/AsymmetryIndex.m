function ai = AsymmetryIndex(theta,TC)

if theta(1) == theta(end)-2*pi
    theta = theta(1:end-1);
    TC = TC(1:end-1);
end


[~,IndMax] = max(TC);
TC = circshift(TC,-(IndMax-round((length(TC)/2))));

% Do an interpolation to smooth the curve
if length(TC)<1e4
    thetaInterp = linspace(0,2*pi,1e4)';
    thetaInterp = thetaInterp';
    TCInterp = interp1(theta',TC,thetaInterp,'spline');
else
    TCInterp = TC;
end

[M,IndMax] = max(TCInterp);
HalfMax = M/2;
[~,IndHalfMaxLeft] = min(abs(TCInterp(1:IndMax)-HalfMax));
[~,IndHalfMaxRight] = min(abs(TCInterp(IndMax:end)-HalfMax));

IndHalfMaxRight = IndHalfMaxRight + IndMax-1;

ai = (thetaInterp(IndMax)-thetaInterp(IndHalfMaxLeft))/...
        (thetaInterp(IndHalfMaxRight)-thetaInterp(IndMax));

