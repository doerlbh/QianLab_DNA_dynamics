% by Baihan Lin, August 2016

function fp = createNatPolymer(fnode, Pc, Pt)
% To create a natural polymer with node nodes

fp = zeros(1,fnode);
fp(2) = 0;

% for no = 3:fnode-1
parfor no = 3:fnode-1
    r = rand();
    if r > Pc*2
        fp(no) = 0;      % not flip, stay trans
    else
        if Pc < r < 2*Pc
            fp(no) = 1;     % flip to cis CCW
        else
            fp(no) = -1;     % flip to cis CW
        end
    end
end

disp('------Simulation-Starts-------');
disp(strcat('createRandPolymer: ', num2str(fp)));
disp('--trials--');
end