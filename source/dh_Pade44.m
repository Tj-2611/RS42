function [h_pade, dh] = dh_Pade44(y,nu,rm,rp,p1,p2,p3,p4,q1,q2,q3,q4)
% rational approximation method by Gatheral and Radoicic (2024)
    h_pade = (p1*y + p2*y.^2 + p3*y.^3 + p4*y.^4) ./ (1 + q1*y + q2*y.^2 + q3*y.^3 + q4*y.^4);
    dh  = 1/2 * (nu*h_pade-rm) .* (nu*h_pade-rp);

end