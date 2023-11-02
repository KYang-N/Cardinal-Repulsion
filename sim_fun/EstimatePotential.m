function Potential = EstimatePotential(Velocity,Location)
% Location is in radian (twice the actual stimulus)
Location = Location/2/pi*180;
Potential = zeros(length(Location),1);
for ii = 2:length(Potential)
Potential(ii) = Potential(ii-1)-Velocity(ii-1)*(Location(ii)-Location(ii-1));
end
Potential = Potential-min(Potential);

