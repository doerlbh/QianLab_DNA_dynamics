% t = 0 ;
% x = 0 ;
% startSpot = 0;
% interv = 1000 ; % considering 1000 samples
% step = 1 ; % lowering step has a number of cycles and then acquire more data
% while ( t <interv )
%     b = sin(t)+5;
%     x = [ x, b ];
%     plot(x) ;
%       if ((t/step)-500 < 0)
%           startSpot = 0;
%       else
%           startSpot = (t/step)-500;
%       end
%       axis([ startSpot, (t/step+50), 0 , 10 ]);
%       grid
%       t = t + step;
%       drawnow;
%       pause(0.003)
% end
  
fig1 = figure;

t = 1 ;
x = 1:500 ;
y = sin(x)+5;
interv = 500 ; % considering 1000 samples
step = 2 ; % lowering step has a number of cycles and then acquire more data
while ( t <interv )
    xt = x(1:t);
    yt = y(1:t);
    plot(xt,yt) ;
      axis([ 0, 500, 0 , 10 ]);
      grid
      t = t + step;
      drawnow;
      pause(0.002)
end
  
saveas(gcf, '/Users/DoerLBH/Dropbox/git/QianLab_DNA_dynamics/test.png','png');
