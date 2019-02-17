function CalculateR

global q m h hhat e epsilon0 u0
clc;
close all;
h = 6.626196e-34;
hhat = h / (2 * pi);
m = 9.109e-31;
e = -1.602e-19;
q = abs(e);
epsilon0 = 8.854187817e-12;
u0 = 4e-7 * pi;  
f = 10e9;
w = 2 * pi * f;
c = 3e8;
lambda = c / f;
k = 2 * pi / lambda;

B = 1e-9;

Rem = c / w;


Pdbm = 30 : -5 : -90;
Pdbm_2E0 = 10 .^ ((Pdbm + 107 - 20 * log10(lambda / 2) + 1.65) / 20.0 - 6);
Bem = Pdbm_2E0 / c;
vEM = Pdbm_2E0 * q / (m * w);
R = Pdbm_2E0 / (w * B);
figure;
plot(Pdbm, log10(Pdbm_2E0), 'b-*');
xlabel('PdBm')
ylabel('log10(E)')
figure;
plot(-Pdbm, log10(R), 'r-*');
figure;
plot(Pdbm, Bem, 'r-*');
L = 0.2;
v0_least = (L * vEM) / (2 * pi * R);
plot(Pdbm, v0_least, 'g-*');

v0 = 1e6;
figure;
plot(Pdbm, log10((L / v0)^2 * sqrt(2) * Pdbm_2E0 * q / m));




