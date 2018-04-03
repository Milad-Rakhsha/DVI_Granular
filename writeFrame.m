function [] = writeFrame(Pos, Vel,currStep)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

global params;
numParticles = size(Pos,1)/2;
headers = {'x', 'y','Vx','Vy','|V|'};
data = zeros(numParticles,5);
filename = strcat('data/state_step.csv.', num2str(currStep));
for i = 1:numParticles
    data(i,1) = Pos(2*i-1);
    data(i,2) = Pos(2*i);
    data(i,3) = Vel(2*i-1);
    data(i,4) = Vel(2*i);
    data(i,5) = norm(Vel(2*i-1:2*i),2);

end

csvwrite_with_headers(filename, data, headers);

end

