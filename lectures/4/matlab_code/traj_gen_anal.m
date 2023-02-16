P1 = [3.4851 2.9752];
P2 = [3.0634 2.9998];
P3 = [1.9537 2.7724];

theta_1 = deg2rad(101.0887);
theta_2 = deg2rad(112.5479);
theta_3 = deg2rad(146.9254);

phi_1 = deg2rad(20);
phi_2 = deg2rad(30);
phi_3 = deg2rad(40);

psi_1 = deg2rad(82.6254);
psi_2 = deg2rad(90.6988);
psi_3 = deg2rad(112.4628);

theta12 = theta_2 - theta_1;
theta13 = theta_3 - theta_1;

phi12 = phi_2 - phi_1;
phi13 = phi_3 - phi_1;

psi12 = psi_2 - psi_1;
psi13 = psi_3 - psi_1;


d12 = P2-P1;
d12 = complex(d12(1),d12(2));

d13 = P3-P1;
d13 = complex(d13(1),d13(2));


a1 = det([d12 exp(theta12*1i)-1; d13 exp(theta13*1i)-1])/det([exp(phi12*1i)-1 exp(theta12*1i)-1; exp(phi13*1i)-1 exp(theta13*1i)-1]);
g1 = det([exp(phi12*1i)-1 d12; exp(phi13*1i)-1 d13])/det([exp(phi12*1i)-1 exp(theta12*1i)-1; exp(phi13*1i)-1 exp(theta13*1i)-1]);

b1 = det([d12 exp(theta12*1i)-1; d13 exp(theta13*1i)-1])/det([exp(psi12*1i)-1 exp(theta12*1i)-1; exp(psi13*1i)-1 exp(theta13*1i)-1]);
f1 = det([exp(psi12*1i)-1 d12; exp(psi13*1i)-1 d13])/det([exp(psi12*1i)-1 exp(theta12*1i)-1; exp(psi13*1i)-1 exp(theta13*1i)-1]);

g1 = 0;
f1 = 0;

z = complex(P1(1),P1(2)) - a1 - g1;
d = complex(P1(1),P1(2)) - z - b1 - f1;
c1 = d + b1 - a1;

a1_m = abs(a1)
g1_m = abs(g1)
b1_m = abs(b1)
f1_m = abs(f1)
z1_m = abs(z)
d1_m = abs(d)
c1_m = abs(c1)


