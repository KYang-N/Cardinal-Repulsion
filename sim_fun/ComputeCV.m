function CV = ComputeCV(Matrix,UnitNum)
% Compute coefficient of variation for all the columns of Matrix
CV = zeros(1,size(Matrix,2));
for kk = 1:length(CV)
    if kk<=8
        ValidPos = UnitNum(:,kk)>=10 & ~isnan(Matrix(:,kk)); % Sessions with excessively few neurons are excluded
        CV(kk) = std(Matrix(ValidPos,kk),0,1)/mean(Matrix(ValidPos,kk),1);
    else
        CV(kk) = std(Matrix(:,kk),0,1,'omitnan')./mean(Matrix(:,kk),1,'omitnan');
    end
end
CV = abs(CV);
CV(CV == 0) = NaN;
end