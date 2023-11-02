% Compute Fourier series
function out = InvFourier(a,b,theta,N)
out = a(1)*ones(1,length(theta));
    for i = 1:N
        out = out+a(i+1)*cos(i*theta)+b(i+1)*sin(i*theta);
    end
end