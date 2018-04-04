function [ B, phi, contact_pairs] = CalcJacobian(NeighborList,Pos )
global params;
nb=params.numParticles;
contact_pairs=[ 0 0];
Jacobian=[];
phi=[];
nc=0;

for i = 1:nb
    
    TOP=1; RIGHT=2; BOTTOM=3; LEFT=4;
    Wall_x=[Pos(2*i-1),params.boxWidth,Pos(2*i-1),0];
    Wall_y=[params.boxHeight,Pos(2*i),0,Pos(2*i)];
    dist_wall=sqrt((Pos(2*i-1)-Wall_x).^2+(Pos(2*i)-Wall_y).^2);
    Normals=[
        0,-1,0,1;
        -1,0,1,0];
    wall_c=dist_wall<params.diam;
    for i_wall=1:4
        if(wall_c(i_wall)==1)
%             fprintf('dist(%d,wall %d)=%1.4f \n',i,i_wall,dist_wall(i_wall)-params.diam/2);
            Dc=zeros(2*(nb+4),1); % Jacobian submatrix
            Dc(2*i-1:2*i,:)=Normals(:,i_wall);
            Dc(2*(nb+i_wall)-1:2*(nb+i_wall),:)=-Normals(:,i_wall);
            Jacobian=[Jacobian Dc];
            contact_pairs=[contact_pairs; [i nb+i_wall]];
            phi=[phi;(dist_wall(i_wall)-params.diam)];
            nc=nc+1;
        end
    end
    
    NL = NeighborList{i};
    m = size(NL,1);
    for neighborIterate = 2:m
        j = NL(neighborIterate);
        dist=sqrt((Pos(2*j-1)-Pos(2*i-1))^2+(Pos(2*j)-Pos(2*i))^2);
        if(dist<params.diam && (sum(contact_pairs(:,1)==j & contact_pairs(:,2)==i)==0))
%             fprintf('dist(%d,%d)=%.4f\n',i,j,dist-params.diam);
            Dc=zeros(2*(nb+4),1); % Jacobian submatrix
            dx=[Pos(2*j-1)-Pos(2*i-1);
                Pos(2*j)-Pos(2*i)];
            n=dx/norm(dx,2);
            v=[ n(2);
                -n(1)];
            Dc(2*i-1:2*i,:)=-n;%[n v];
            Dc(2*j-1:2*j,:)=n;%[n v];
            Jacobian=[Jacobian Dc];
            contact_pairs=[contact_pairs; [i j]];
            phi=[phi;(n'*dx-params.diam)];
            nc=nc+1;
        end
    end
    
end
contact_pairs=contact_pairs(2:end,1:2);
B=Jacobian;

end