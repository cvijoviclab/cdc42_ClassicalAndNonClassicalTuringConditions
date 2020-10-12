function [uStar,vStar, VStar] = SymbolicSS(c1,c_1,c2,a,V0_init,cmax)
%SombolicSS Uses a symbolic approach to finding the steady states
%   Detailed explanation goes here

syms u
syms v

m= V0_init/a;

eqns = [c2*v - u + u*u*v == 0, c1*(V0_init-a*(u+v)).*(cmax-u-v) - c_1*v == 0];
sol = solve(eqns,[u v],'ReturnConditions',true);
U = vpa(sol.u);
V = vpa(sol.v);
k = length(U);
removeind = [];

for i = 1:k
    if ~isreal(U(i)) || U(i) < 0 || V(i) < 0 || U(i) + V(i) > min(cmax,m)
        removeind = [removeind, i];
    end
end
U(removeind) = [];
V(removeind) = [];

% Gather all steady states
uStar = U;
vStar = V;
VStar = (V0_init.*ones(size(uStar)))-(a.*(uStar+vStar));


end

