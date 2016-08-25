% by Baihan Lin, August 2016

function [fACorT] = AutocorEnd(ftP)
% calculate autocorrelation for t

n = size(ftP,2);
fACorT = zeros(1, n);
Cst = ftP(:,1);

for c = 1:n
    Cfin = ftP(:,c);
    fACorT(c) = dot(Cst, Cfin);
end

end
