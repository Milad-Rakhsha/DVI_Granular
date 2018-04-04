function [newPos,newVel] = integrator(Pos,Vel)
global params;
Pos=[Pos;zeros(8,1)];
Vel=[Vel;zeros(8,1)];


[ NeighborList ] = GetNeighborList( Pos, params.diam, params.boxWidth, params.boxHeight,params.numParticles );
[ B, phi, contact_pairs ] = CalcJacobian(NeighborList,Pos);
nb=params.numParticles+4;
nc=size(contact_pairs,1);
MInv=(eye(2*nb)*1/params.particleMass);
F_ext=repmat([0 ; params.particleMass*params.g],params.numParticles,1);
F_ext=[F_ext; zeros(8,1)];
F_contact=F_ext*0;
if(size(B,2)~=0)
N=B'*MInv*B;
p=phi/params.dt+B'*Vel+params.dt*B'*MInv*F_ext;
A=-eye(nc,nc);
b=zeros(nc,1);
options = optimoptions('quadprog','Algorithm','interior-point-convex',...
    'Display','off');
% f = quadprog(N,p,A,b,[],[],[],[],[],options);
f = IPM(N,p,params);

F_contact=B*f;
end

newVel=MInv*params.dt*(F_contact+F_ext)+Vel;

newPos=Pos+params.dt*newVel;
newPos=newPos(1:params.numParticles*2);
newVel=newVel(1:params.numParticles*2);

