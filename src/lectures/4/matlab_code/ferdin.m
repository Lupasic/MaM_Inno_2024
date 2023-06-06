theta_1 = deg2rad(173.9);
theta_2 = deg2rad(83.9);
theta_3 = deg2rad(141.1958);
psi_1 = deg2rad(7.6);
psi_2 = deg2rad(72.5);
psi_3 = deg2rad(60.1145);

a = 1;

A = [cos(psi_1) -cos(theta_1) 1 ;
    cos(psi_2) -cos(theta_2) 1;
    cos(psi_3) -cos(theta_3) 1];
B = [cos(psi_1 - theta_1); cos(psi_2 - theta_2); cos(psi_3 - theta_3)];

X = linsolve(A,B)

syms b c d;
eqns = [X(1) == d/a; X(2) == d/b; X(3) == (d^2 + b^2 + a^2 - c^2)/(2*a*b)];
res = solve(eqns);
abs(double(res.b))
abs(double(res.c))
abs(double(res.d))



