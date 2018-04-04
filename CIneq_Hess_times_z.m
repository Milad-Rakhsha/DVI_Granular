function [out] = CIneq_Hess_times_z(x,z)
H_Ineq=zeros(length(x));
for i=1:length(z)
H_Ineq=H_Ineq+ zeros(length(x))*z(i);
end
out =H_Ineq;

end