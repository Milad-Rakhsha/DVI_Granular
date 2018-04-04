function [Pos,Vel,Acc,flag] = initParticleSystem(params)
%initParticleSystem 
%   This function initialize an array of particles.

disp('Initializing Particle System');

numParticles=params.numParticles;
numParticlesPerRow = params.numParticlesPerRow;
numParticlesPerCol = params.numParticlesPerCol;
dr = params.initialSeparation;

Pos = zeros(2*numParticles,1);
Vel = zeros(2*numParticles,1);
Acc = zeros(2*numParticles,1);
counter = 1;

% Place Particles
for iy = 0: numParticlesPerCol-1
    for ix = 0:numParticlesPerRow-1 
        Pos(counter) = ix*dr+ dr/2+4*dr;
        Pos(counter+1) = iy*dr*1.2+ dr/2 ;
        counter = counter + 2; 
    end
end
