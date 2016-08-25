% by Baihan Lin, August 2016

function [fPnew, ftP, fHTd, ffin] = fasttwistLoopRand(fpath, ft, fp, fa, fL, fangle, fHc, fHt)
% To twist randomly till formed a loop

disp(strcat('T-',num2str(ft),'-------------'));

fPnew = fp;
ffin = 1;

ftP = finConfig(fPnew, fL, fangle);

% stair = length(fp)-1;
%
% fig = figure;
% fv = visV(buildV(fPnew, fL, fangle));
% xt = fv(1,:);
% yt = fv(2,:);
% zt = fv(3,:);
% plot3(xt(1:stair),yt(1:stair),zt(1:stair));
% xlmin = -length(fp)*fL;
% xlmax = length(fp)*fL;
% ylmin = -length(fp)*fL;
% ylmax = length(fp)*fL;
% zlmin = -length(fp)*fL;
% zlmax = length(fp)*fL;
% axis([ xlmin, xlmax, ylmin, ylmax, zlmin, zlmax]);
% grid;
% xlabel 'x';
% ylabel 'y';
% zlabel 'z';
% title(strcat('3D-Simulation-of-Node-',num2str(length(fp)),'-Trial-',num2str(ft),'-rigid-polymer-dynamics'));
% axis([ xlmin, xlmax, ylmin, ylmax, zlmin, zlmax]);
%
% disp(strcat('Debug ',num2str(ffin),': ', num2str(fPnew)));

while HTdist(fPnew, fL, fangle) > fa
    
    no =  randsample(2:length(fp)-1,1);
    
    Pno = randFlip(fPnew(no));
    Phypo1 = [fPnew(1:no-1), Pno(1), fPnew(no+1:end)];    % hypothetical change
    Phypo2 = [fPnew(1:no-1), Pno(2), fPnew(no+1:end)];    % hypothetical change
    
    Pch1 = pE(Phypo1, fHc, fHt)/(pE(Phypo1, fHc, fHt)+pE(Phypo2, fHc, fHt)+pE(fPnew, fHc, fHt));
    Pch2 = pE(Phypo2, fHc, fHt)/(pE(Phypo1, fHc, fHt)+pE(Phypo2, fHc, fHt)+pE(fPnew, fHc, fHt));
    rtt = rand();
    if rtt < Pch1
        fPnew = Phypo1;      % change state
    else
        if rtt < Pch1+Pch2
            fPnew = Phypo2;
        end
    end
    if HTdist(fPnew, fL, fangle) < fa
        break;
    end
    
    ftP = [ftP finConfig(fPnew, fL, fangle)];
    ffin = ffin+1;
    %
    %     stair = length(fp)-1;
    %
    %     fv = visV(buildV(fPnew, fL, fangle));
    %     xt = fv(1,:);
    %     yt = fv(2,:);
    %     zt = fv(3,:);
    %
    %     plot3(xt(1:stair), yt(1:stair), zt(1:stair));
    %     grid;
    %     xlmin = -length(fp)*fL;
    %     xlmax = length(fp)*fL;
    %     ylmin = -length(fp)*fL;
    %     ylmax = length(fp)*fL;
    %     zlmin = -length(fp)*fL;
    %     zlmax = length(fp)*fL;
    %     axis([ xlmin, xlmax, ylmin, ylmax, zlmin, zlmax]);
    %     xlabel 'x';
    %     ylabel 'y';
    %     zlabel 'z';
    %     title(strcat('3D-Simulation-of-Node-',num2str(length(fp)),'-Trial-',num2str(ft),'-rigid-polymer-dynamics'));
    %     xc = xlim;
    %     xl = xc(1)*0.2+xc(2)*0.8;
    %     yc = ylim;
    %     yl1 = yc(1)*0.14+yc(2)*0.86;
    %     yl2 = yc(1)*0.26+yc(2)*0.74;
    %     zc = zlim;
    %     zl = zc(1)*0.2+zc(2)*0.8;
    % %
    %     text(xl,yl1,zl, strcat('ffin=',num2str(ffin)),'Color','red','FontSize',12);
    %     text(xl,yl2,zl, strcat('HTdist=',num2str(HTdist(fPnew, fL, fangle))),'Color','red','FontSize',12);
    %     drawnow;
    
    %       pause(0.6)
    
    %     disp(strcat('Debug ',num2str(ffin),': ', num2str(fPnew)));
    
end

fHTd = HTdist(fPnew, fL, fangle);

% filename = strcat(fpath, 'N',num2str(length(fp)),'-T',num2str(ft),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
% parsaveas(gcf, filename,'png');
% close gcf;

filename = strcat(fpath, 'Loop-ftP-N',num2str(length(fp)),'-f',num2str(ffin),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.txt');
parsave(filename, ftP, '-ascii');

disp(strcat('final state: ', num2str(fPnew)));
disp(strcat('finish time: ', num2str(ffin)));
disp(strcat('finish dist: ', num2str(HTdist(fPnew, fL, fangle))));

end