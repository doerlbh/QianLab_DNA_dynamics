% by Baihan Lin, August 2016

function fm = buildV(fp, fL, fangle)
% To build a vector set (matrix) based on given polymer states

ftrans = ones(1, length(fp)-2);
for count = 2:length(fp)-2
    ftrans(count) = -ftrans(count-1);
end
fm2D = build2DV([0, ftrans, 0], fL, fangle);

fx = fm2D(1,:);
fy = fm2D(2,:);
fz = zeros(1, length(fp)-1);
% length(fx)
% length(fy)
% length(fz)

fm = [fx;fy;fz];

for no = 2:length(fp)-1
    if fp(no) == 0
    else
        ax = fm(:,no-1)/norm(fm(:,no-1));
        u = ax(1);
        v = ax(2);
        w = ax(3);
        
        if fp(no) == 1
            ang = 2*pi/3; % CCW
        else
            ang = 4*pi/3; % CW
        end
        rot = [u^2+(1-u^2)*cos(ang), u*v*(1-cos(ang))-w*sin(ang), u*w*(1-cos(ang))+v*sin(ang);
            u*v*(1-cos(ang))+w*sin(ang), v^2+(1-v^2)*cos(ang), v*w*(1-cos(ang))-u*sin(ang);
            u*w*(1-cos(ang))-v*sin(ang), w*v*(1-cos(ang))+u*sin(ang), w^2+(1-w^2)*cos(ang)];
        fm(:,no:end) = rot*fm(:,no:end);
    end
end

end



