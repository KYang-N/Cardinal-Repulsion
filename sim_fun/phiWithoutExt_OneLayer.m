function out = phiWithoutExt_OneLayer(W,s,IEc,qm)
Input = W*s+IEc;
out = qm(Input);
end