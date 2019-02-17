
function RingMagneticField
clc;
close all;
format long;
global q m h hhat e epsilon0 u0
h = 6.626196e-34;
hhat = h / (2 * pi);
m = 9.109e-31;
e = -1.602e-19;
q = abs(e);
epsilon0 = 8.854187817e-12;
u0 = 4e-7 * pi;  
f = 12e6;
w = 2 * pi * f;
c = 3e8;
lambda = c / f;
k = 2 * pi / lambda;

L = 0.25;
v0 = 1e6;
simcount = 2000;
tmin = 0;
tmax = L / v0;
t = linspace(tmin, tmax, simcount);
dt = (tmax - tmin) / length(t);
E0 = 50;
vEM = q * E0 / (m * w); 
xEM = q * E0 / (m * w^2);

AXES_COUNT = 3;
global isConvergingMagnetic
isConvergingMagnetic = 0;
if isConvergingMagnetic == 0
    B = 0.03e-3;
else
    B = 0.3e-3;
end
vx = v0 * ones(1, simcount);
vy = zeros(1, simcount);
vz = zeros(1, simcount);
xy = zeros(1, simcount);
xx = zeros(1, simcount);
xz = zeros(1, simcount);

for i = 2: simcount - 1
    xx(i) = xx(i - 1) + dt * (vx(i) + vx(i + 1)) / 2.0;
    xy(i) = xy(i - 1) + dt * (vy(i) + vy(i + 1)) / 2.0;
    xz(i) = xz(i - 1) + dt * (vz(i) + vz(i + 1)) / 2.0;
end

xx(simcount) = xx(simcount - 1);
xy(simcount) = xy(simcount - 1);
xz(simcount) = xz(simcount - 1);

global xVpp
xVpp = max(xy);
rmax = (1.5 * xVpp) / 2;
theta = 0 : pi/100 : 2 * pi;
By = rmax * cos(theta) + xVpp / 2.0; 
Bz = rmax * sin(theta);

figure;
scatter3(xx, xy, xz);
hold on;
global startBx endBx stepBx Barea
startBx = 0.105;
endBx = 0.18;
stepBx = 0.005;
Barea = startBx : stepBx : endBx;

for bx = startBx : stepBx : endBx
    Bx = ones(length(theta), 1) * bx;
    scatter3(Bx, By, Bz, 'r')
end
hold off;

v_xyz = [vx' vy' vz'];

B_xyz = [zeros(simcount, 1), zeros(simcount, 1), [zeros(simcount * startBx / L , 1); ones(round(simcount * (endBx - startBx) / L), 1) * B; zeros(round(simcount * (L - endBx) / L), 1)]];

trajectory_xyz = zeros(simcount, AXES_COUNT);

for i = 2 : simcount * startBx / L    
    i
    E_acc_xyz = [0, q * E0 / m * cos(w * t(i) + 0.5 * pi), 0];
    v_xyz(i + 1, 1:end) = v_xyz(i, 1:end) + (E_acc_xyz) * dt;
    trajectory_xyz(i, 1:end) = trajectory_xyz(i - 1, 1:end) + 0.5 * dt * (v_xyz(i, 1:end) + v_xyz(i + 1, 1:end)); 
end

for i = simcount * startBx / L + 1 : simcount * endBx / L;
    i;
    last_xyz = trajectory_xyz(i - 1, 1:end);
    if isConvergingMagnetic == 0
        if last_xyz(2) > 1e-12
            B_xyz(i, 1:end) = B_xyz(i, 1:end);
        else
            B_xyz(i, 1:end) = [0 0 0];
        end
    else
        if last_xyz(2) < xVpp * 0.5
            B_xyz(i, 1:end) = B_xyz(i, 1:end);
        else
            B_xyz(i, 1:end) = -B_xyz(i, 1:end);
        end
    end
    
    E_acc_xyz = [0 q * E0 / m * cos(w * t(i)) 0];
    if getValueAbs(B_xyz) == 0
        lorentz_acc_xyz = [0 0 0];
    else
        lorentz_acc_xyz = getLorentz(v_xyz, B_xyz, i) / m;
    end
    next_v = v_xyz(i, 1:end) + (lorentz_acc_xyz + E_acc_xyz) * dt;
    % v_xyz(i + 1, 1:end) = next_v * getValueAbs(v_xyz(i, 1:end)) / getValueAbs(next_v);
    v_xyz(i + 1, 1:end) = next_v;
    next_xyz = trajectory_xyz(i - 1, 1:end) + 0.5 * dt * (v_xyz(i, 1:end) + v_xyz(i + 1, 1:end));
    trajectory_xyz(i, 1:end) = next_xyz;
end

for i = simcount * endBx / L + 1 : simcount - 1;
    i
    E_acc_xyz = [0 q * E0 / m * cos(w * t(i)) 0];
    next_v = v_xyz(i, 1:end) + E_acc_xyz * dt;
    v_xyz(i + 1, 1:end) = next_v;
    next_xyz = trajectory_xyz(i - 1, 1:end) + 0.5 * dt * (v_xyz(i, 1:end) + v_xyz(i + 1, 1:end));
    trajectory_xyz(i, 1:end) = next_xyz;
end

trajectory_xyz(simcount, 1:end) = trajectory_xyz(simcount - 1, 1:end);
showTrajectoryXY(trajectory_xyz);
showVxyzAndTime(v_xyz, t);

function vet = getTimevalue(values, t)
vet = values(t, 1:end);

function valueabs = getValueAbs(value)
valueabs = sqrt(sum(value .* value));

function lorentz = getLorentz(vv, BB, index)
global e
if nargin == 2
    lorentz = e * cross(vv, BB);
end
if nargin == 3
    lorentz = e * cross(getTimevalue(vv, index), getTimevalue(BB, index));
end

function showTrajectoryXYZ(trajectory_xyz)
figure;
scatter3(trajectory_xyz(1:end, 1), trajectory_xyz(1:end, 2), trajectory_xyz(1:end, 3));

function showTrajectoryXY(trajectory_xyz)
figure;
scatter(trajectory_xyz(1:end, 1), trajectory_xyz(1:end, 2));
hold on
global startBx endBx stepBx xVpp isConvergingMagnetic
for x = startBx + stepBx / 2 : stepBx : endBx
    for y = -xVpp * 3.0 : stepBx / 2: xVpp * 3.0
        vpp = y + xVpp * 0.5;
        if isConvergingMagnetic == 0
            if vpp >= xVpp * 0.5
                scatter(x, vpp, 'r')
            else
                scatter(x, vpp, 'r*')
            end
        else
            if vpp < xVpp * 0.5
                scatter(x, vpp, 'r')
            else
                scatter(x, vpp, 'r*')
            end
        end
    end
end
plot(linspace(0, 0.2, 50), ones(50, 1) * xVpp / 2.0, 'b');
hold off

function showVxyzAndTime(vxyz, t)
v = sqrt(vxyz(1:end, 1).^2 + vxyz(1:end, 2).^2 + vxyz(1:end, 3).^2);
figure;
plot(t, v);
