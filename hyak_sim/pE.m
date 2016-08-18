
function fE = pE(fp, fHc, fHt)
% To calculate energy of a certain polymer state

fE = 0;

% for no = 3:length(fp)-1
parfor no = 3:length(fp)-1
    if abs(fp(no)) == 1
        fE = fE + fHc;      % cis
    else
        fE = fE + fHt;     % trans
    end
end

end

