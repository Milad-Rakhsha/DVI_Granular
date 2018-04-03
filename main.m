clear
clc
mkdir data
delete data/state_step.csv.*
global params

params.numSteps = 1000; % Number of frames
params.dt = 1e-2; % Time step
params.diam = 0.1; % particle diameter
params.initialSeparation = params.diam;
params.rho0 = 1000;
params.particleMass = pi*(params.diam ^2/4) *params.rho0;
params.g = -0.01; % Gravity acceleration
params.boxWidth = 20 * params.diam;
params.boxHeight = 20 * params.diam ;
params.numParticlesPerRow =11;
params.numParticlesPerCol =11;
params.numParticles = params.numParticlesPerRow*(params.numParticlesPerCol); %Filename

currStep = 0;
[Pos, Vel, Acc] = initParticleSystem(params);
drawFrame(Pos, params.numParticles, params.boxWidth, params.boxHeight, params.diam);

counter=0;
fprintf('Simulation loop\n')
for i=1:1e6
    [newPos,newVel] = integrator(Pos,Vel);
    Pos=newPos;
    Vel=newVel;
    if(mod(i,50)==0)
        fprintf('---------Frame %d\n',i);
        writeFrame(newPos, newVel,counter);
        counter=counter+1;
    end
    
end
