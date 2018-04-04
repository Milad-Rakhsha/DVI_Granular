clear
clc
mkdir data
delete data/state_step.csv.*
global params

params.numSteps = 1e6; % Number of frames
params.dt = 1e-3; % Time step
params.diam = 0.05; % particle diameter
params.initialSeparation = params.diam;
params.rho0 = 1000;
params.particleMass = pi*(params.diam ^2/4) *params.rho0;
params.g = -0.01; % Gravity acceleration
params.boxWidth = 40 * params.diam;
params.boxHeight = 40 * params.diam ;
params.numParticlesPerRow =21;
params.numParticlesPerCol =21;
params.numParticles = params.numParticlesPerRow*(params.numParticlesPerCol); %Filename
% Solver settings
params.maxNewtownIter=20;
params.maxtotalIter=8;
params.epsilon=1e-10;
params.taw=0.995;

currStep = 0;
[Pos, Vel, Acc] = initParticleSystem(params);
% drawFrame(Pos, params.numParticles, params.boxWidth, params.boxHeight, params.diam);

counter=0;
fprintf('Simulation loop\n')
for i=1:params.numSteps
    [newPos,newVel] = integrator(Pos,Vel);
    Pos=newPos;
    Vel=newVel;
    if(mod(i,100)==0)
        fprintf('---------Frame %d\n',i);
        writeFrame(newPos, newVel,counter);
        counter=counter+1;
    end
    
end
