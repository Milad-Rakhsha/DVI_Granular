function [xs] = IPM(H,f,params)

% This function solves an optimization problem for :
% min 1/2*x'*H*x+x'*f s.t. CIneq>=0;
% This is posed as a linear system of the form (Page 570 Nocedel and Wright)
% [ Hess_xx(L)+A_Ineq^T*(S^-1*Z)*A_Ineq ] [Dx] = -[Grad(f)-A_Ineq'(x)*z + A_Ineq'*(inv(S)*Z)*(CIneq-mu*inv(Z)*e)]
% with         Z*ds+S*dz=-(S*z-mu*e);
% AND          dz=-(S^-1*Z)*(CIneq-mu*Z^-1*e+A_Ineq*dx)
% Sigma=S^-1*Z
% x is the input variable (n*1)
% where Hess_xx(L) is the Hessian of the Lagrangian function f-CIneq^T(x)*z (n*n);
% Grad(f) is the gradient of the object function (n*1).
% A_Ineq(x) is the jacobian of the Ineq
% constraints w.r.t. x and is respectively n_Ineq*n; y and z are the lagrange multiplier associated with
% the CIneq,so they have the same size as the number of Ineq constraints (n_Ineq*1).
% Y, Z are basically y and z that are placed on the diagonal of a matrix :  (n_Ineq*n_Ineq)
% s is the slack variable to convert the inequality constraints into the
% standard format (n_Ineq*1). S is the diagonal matrix build from s (n_Ineq*n_Ineq)
n=length(f);
x=zeros(n,1);
n_Ineq=length(CIneq(x));
mu=10;
s=CIneq(x)+0.01;
z=mu./s;
e=ones(n_Ineq,1);

% fprintf('n=%d, n_Ineq=%d\n',n,n_Ineq);
% fprintf('\n\t\t#Newton     Res      alpha_s   alpha_z     mu\n');
% fprintf('\t\t-------------------------------------------------------------------------------------\n');
mu=s'*z/n_Ineq;
for Newton=1:params.maxtotalIter
    alpha_s=0.95; alpha_z=0.95;
    for k=1:Newton
        H_Lag=H-CIneq_Hess_times_z(x,z);
        A_Ineq=CIneq_grad(x);
        S=diag(s);
        Z=diag(z);
        A_Matrix=H_Lag+A_Ineq'*(diag(z./s))*A_Ineq;
        b=-((H*x+f)-A_Ineq'*z+...
           A_Ineq'*(diag(z./s))*(CIneq(x)-mu*diag(1./z)*e));
        
        dx=A_Matrix\b;
        dz=diag(z./s)*(-(CIneq(x)-mu*diag(1./z)*e)-A_Ineq*dx);
        ds=diag(1./z)*(-(S*z-mu*e)-S*dz);

        [alpha_s, alpha_z] = FindAlphas(x,s,z,dx,ds,dz,params.taw);
        x=x+alpha_s*dx;
        z=z+alpha_z*dz;
        s=s+alpha_s*ds;
        Error=norm(b,inf);
        if(Error<params.epsilon)
            break;
        end
    end
    mu=mu*0.1;
end
%     fprintf('Iter %d\t\t %d\t %2.2e, %2.2e, %2.2e, %2.2e\n',Newton,k,norm(b,inf),alpha_s,alpha_z,mu);
    xs=x;

end


