function [results] = padeLaplace(time,decay,irf)

%I AM NOT COMPLETELY SURE TAYLOR SERIES CALCULATIONS ARE CORRECT

p0 = 1;
nTay = 10;
x = zeros(nTay,1); % x is the list of coefficients for the taylor series

%determine the taylor series coefficients by the trapezoid rule
for i = 0:nTay
    timeTerm = (-1*time).^i;
    decayExpTerm = decay.*exp(-p0*time);    
    x(i+1,1) = 1/factorial(i)*trapz(time,timeTerm.*decayExpTerm);
end

%results shows the poles and residues for each "level" of pade-laplace
results = struct('L_M',"",'n',0,'d',0,'poles',0,'residues',0);

%calculate the relevant pade coefficients (this is brute force not
%recursive). note that n0 and d0 are already defined
%[0/1]
results(1,1).L_M = "0_1";
%using the identity that n0 = x0
results(1,1).n = x(1,1);
d1 = -1*x(2,1)/x(1,1);
results(1,1).d = [1 d1];

%[1/2], answers array will be in the order [n1 d1 d2]
syms n1 d1 d2
eqn1 = x(2,1) + x(1,1)*d1 == n1;
eqn2 = x(3,1) + x(2,1)*d1 + x(1,1)*d2 == 0;
eqn3 = x(4,1) + x(3,1)*d1 + x(2,1)*d2 == 0;
[A,B] = equationsToMatrix([eqn1, eqn2, eqn3], [n1, d1, d2]);
val = linsolve(A,B);
results(2,1).n = [x(1,1) val(1)];
results(2,1).d = [1 val(2) val(3)]; %still symbolic, hopefully this is OK

%[2/3], answers array will be in the order [n1 n2 d1 d2 d3]

%[4/5], answers array will be in the order [n1 n2 n3 d1 d2 d3 d4]

%[5/6], answers array will be in the order [n1 n2 n3 n4 d1 d2 d3 d4 d5]

end

