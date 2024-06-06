function DecodedOrientation = PVDecoder(PreferredFeature,FiringRate)
% PreferredFecture is a row vector and FiringRate is a column vector
if max(FiringRate)<0.1
    DecodedOrientation = nan(1,size(FiringRate,2));
else
    DecodedOrientation = angle(exp(1i*PreferredFeature)*FiringRate/sum(FiringRate));
    if DecodedOrientation < 0
        DecodedOrientation = DecodedOrientation + 2*pi;
    end
end