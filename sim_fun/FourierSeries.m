% Use Fourier series to fit the velocity and coefficient of diffusion
function [a,b] = FourierSeries(data,NFourier,ThetaLoc)
dtheta = diff(ThetaLoc);
dtheta(dtheta<0) = dtheta(dtheta<0) + 2*pi;
ThetaLoc = ThetaLoc(1:end-1);

data = data(:,1:end-1);
cosn = zeros(NFourier,length(ThetaLoc));
sinn = cosn;
for i = 1:NFourier
    cosn(i,:) = cos(i*ThetaLoc);
    sinn(i,:) = sin(i*ThetaLoc);
end
a = zeros(1,NFourier+1);
b = zeros(1,NFourier+1);
a(1) = 1/(2*pi)*data*dtheta';
b(1) = 0;
for i = 1:NFourier
    a(i+1) = 1/pi*cosn(i,:)*(data.*dtheta)';
    b(i+1) = 1/pi*sinn(i,:)*(data.*dtheta)';
end
end