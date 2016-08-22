% by Baihan Lin, August 2016

function fC = finConfig(fp, fL, fangle)
% To calculate the final configuration of end to end of a polymer

fm = buildV(fp, fL, fangle);
fC = sum(fm,2);

end