function PF = ComputePF(TC,SampleInput)
% TC should be a NxM matrix, where N denotes the number of neurons and M
% denotes the number of input orientations. SampleInput is a row vector.

InputInterpNum = 1e3;
if length(SampleInput) < InputInterpNum
    InputInterp = linspace(0,2*pi,InputInterpNum);
    TCInterp = zeros(size(TC,1),InputInterpNum);
    for i = 1:size(TC,1)
        TCInterp(i,:) = interp1(SampleInput,TC(i,:),InputInterp,'spline');
    end
end
[~,Ind] = max(TCInterp,[],2);
PF = InputInterp(Ind);
PF(1) = 0;