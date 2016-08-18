
function fp = createRandPolymer(fnode)
% To create a random polymer with node nodes

fp = zeros(1,fnode);
fp(2) = 0;

% for no = 3:fnode-1
parfor no = 3:fnode-1
    r = rand();
    if r < 1/3
        fp(no) = 0;      % not flip, stay trans
    else
        if 1/3 < r < 2/3
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