function value = q(u,v,V0_init,a,cmax,c1,c_1)
 value = ( c1 * (cmax - (u+v)) * (V0_init - (a * (u+v))) ) - (c_1 * v);
end