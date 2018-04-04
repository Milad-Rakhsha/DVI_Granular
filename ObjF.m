function [out] = ObjF(x,H,f)
out=1/2*x'*H*x+x'*f;
end