function DeltaV = HCW_DeltaV(p0,pt,v0,vt,t,omega)

%% Phi matricies 

phi_11 = [4-3*cos(omega*t) 0 0;
          6*sin(omega*t)-6*omega*t 1 0;
          0 0 cos(omega*t)];
   
phi_21 = [3*omega*sin(omega*t) 0 0;
          6*omega*cos(omega*t)-6*omega 0 0;
          0 0 -omega*sin(omega*t)];
      
phi_22 = [cos(omega*t) 2*sin(omega*t) 0;
          -2*sin(omega*t) 4*cos(omega*t)-3 0;
          0 0 cos(omega*t)];

phi_12_inv = [ (omega*(4*sin(omega*t) - 3*omega*t))/(4*cos(omega*t)^2 - 8*cos(omega*t) + 4*sin(omega*t)^2 - 3*omega*t*sin(omega*t) + 4), (2*omega*(cos(omega*t) - 1))/(4*cos(omega*t)^2 - 8*cos(omega*t) + 4*sin(omega*t)^2 - 3*omega*t*sin(omega*t) + 4), 0;
 -(2*omega*(cos(omega*t) - 1))/(4*cos(omega*t)^2 - 8*cos(omega*t) + 4*sin(omega*t)^2 - 3*omega*t*sin(omega*t) + 4), (omega*sin(omega*t))/(4*cos(omega*t)^2 - 8*cos(omega*t) + 4*sin(omega*t)^2 - 3*omega*t*sin(omega*t) + 4), 0;
0, 0, omega/sin(omega*t)];


DeltaV1 = phi_12_inv*(pt-phi_11*p0)-v0;
DeltaV2 = vt - (phi_21*p0 + phi_22*phi_12_inv*(pt-phi_11*p0));

DeltaV = abs(DeltaV1) + abs(DeltaV2);
end



