function [fACorT] = AutocorEnd(ftP)
% calculate autocorrelation for t

n = length(ftP);
fACorT = zeros(1, n);
Cst = ftP(1);
Re0 = dot(Cst, Cst);

for c = 1:n
    Cfin = ftP(n);
    fACorT(c) = dot(Cst, Cfin)/Re0;
end

end
