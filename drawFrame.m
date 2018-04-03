function [] = drawFrame(Pos, numParticles, boxWidth, boxHeight, radius)
%drawFrame Plot particles at a specific time step
hold on;
% axis([-0.2 boxWidth+0.5 -0.2 boxHeight+0.5]);

for i = 1:numParticles
    x = Pos(2*i-1);
    y = Pos(2*i);
    scatter(x,y,'MarkerEdgeColor',[0 .5 .5],...
        'MarkerFaceColor',[0 .7 .7],...
        'LineWidth',1.5);
end


hold off;


end

