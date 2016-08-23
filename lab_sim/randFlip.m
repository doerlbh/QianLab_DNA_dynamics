% by Baihan Lin, August 2016

function fPno = randFlip(fn)
% To randomly flip

if fn == 0
    fPno = [1, -1];
else
    fPno = [0, -fn];
end

end

