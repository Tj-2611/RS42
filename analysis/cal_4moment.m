function [Ex, Var, Skew, Kurt] = cal_4moment(params, model, du, IR, t)

E0 = 1;
n = 1e2;
Edu = charfunv2(model,params,n,du*(-1i),t);
E2du = charfunv2(model,params,n,2*du*(-1i),t);
E3du = charfunv2(model,params,n,3*du*(-1i),t);
E4du = charfunv2(model,params,n,4*du*(-1i),t);

E_du = charfunv2(model,params,n,-du*(-1i),t);
E_2du = charfunv2(model,params,n,-2*du*(-1i),t);
E_3du = charfunv2(model,params,n,-3*du*(-1i),t);
E_4du = charfunv2(model,params,n,-4*du*(-1i),t);

E0 = log(real(E0));
Edu = log(real(Edu)) + IR*t*du;
E2du = log(real(E2du)) + IR*t*2*du;
E3du = log(real(E3du)) + IR*t*3*du;
E4du = log(real(E4du)) + IR*t*4*du;
E_du = log(real(E_du)) + IR*t*(-1)*du;
E_2du = log(real(E_2du)) + IR*t*(-2)*du;
E_3du = log(real(E_3du)) + IR*t*(-3)*du;
E_4du = log(real(E_4du)) + IR*t*(-4)*du;



d1E0 = (Edu-E_du)/(2*du);
d1Edu = (E2du-E0)/(2*du);
d1E_du = -(E_2du -E0)/(2*du);
d1E2du = (E3du-Edu)/(2*du);
d1E_2du = -(E_3du-E_du)/(2*du);
d1E3du = (E4du-E2du)/(2*du);
d1E_3du = -(E_4du-E_2du)/(2*du);

d2E0 = (d1Edu - d1E_du)/(2*du);
d2Edu = (d1E2du - d1E0)/(2*du);
d2E_du = -(d1E_2du - d1E0)/(2*du);
d2E2du = (d1E3du - d1Edu)/(2*du);
d2E_2du = -(d1E_3du - d1E_du)/(2*du);

d3E0 = (d2Edu-d2E_du)/(2*du);
d3Edu = (d2E2du-d2E0)/(2*du);
d3E_du = -(d2E_2du-d2E0)/(2*du);


d4E0 = (d3Edu-d3E_du)/(2*du);

Ex = d1E0;
Var = d2E0;
Skew = d3E0/((d2E0)^(3/2));
Kurt = d4E0/(d2E0)^2;


end


