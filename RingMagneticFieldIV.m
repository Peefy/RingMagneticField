
function RingMagneticFieldIV

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
f = 2e6;
w = 2 * pi * f;
c = 3e8;
lambda = c / f;
v0 = 1e6;
simcount = 2500;
B = 0.3e-1;
B_xyz = [ones(simcount, 1) * 0, ones(simcount, 1) * 0, ones(simcount, 1) * B];
tmin = 0;
r0 = m * v0 / (q * B);
tmax = r0 * 2 * pi / v0 * 2;
t = linspace(tmin, tmax, simcount);
dt = (tmax - tmin) / length(t);
AXES_COUNT = 3;

vx = zeros(1, simcount);
vy = zeros(1, simcount);
vz = zeros(1, simcount);
xx = zeros(1, simcount);
xy = zeros(1, simcount);
xz = zeros(1, simcount);

vy(1) = v0;
vy(2) = v0;
xx(1) = r0;
xx(2) = r0;

v_xyz = [vx' vy' vz'];

trajectory_xyz = [xx' xy' xz'];


for i = 2 : simcount - 1;
    i
    last_xyz = trajectory_xyz(i - 1, 1:end);
    last_r = sqrt(last_xyz(1)^2 + last_xyz(2)^2);
    last_phi = atan(last_xyz(2) / last_xyz(1));
    Bxyz_now = B_xyz(i, 1:end);
    %Bxyz_now(1) = 8 * B * cos(last_phi) / ((last_r + r0) / r0)^3;
    %Bxyz_now(2) = 8 * B * sin(last_phi) / ((last_r + r0) / r0)^3;
    vxyz_now = v_xyz(i, 1:end);
    if getValueAbs(B_xyz) == 0 
        lorentz_acc_xyz = [0 0 0];
    else
        lorentz_acc_xyz = getLorentz(vxyz_now, Bxyz_now) / m;
    end 
    next_v = v_xyz(i, 1:end) + (lorentz_acc_xyz) * dt;
    v_xyz(i + 1, 1:end) = next_v;
    next_xyz = trajectory_xyz(i - 1, 1:end) + 0.5 * dt * (v_xyz(i, 1:end) + v_xyz(i + 1, 1:end));
    trajectory_xyz(i, 1:end) = next_xyz;
end

trajectory_xyz(simcount, 1:end) = trajectory_xyz(simcount - 1, 1:end);
showTrajectoryXYZWithoutB(trajectory_xyz);

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
    v = getTimevalue(vv, index);
    b = getTimevalue(BB, index);
    lorentz = e * cross(v, b);
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

function showTrajectoryXYWithoutB(trajectory_xyz)
figure;
scatter(trajectory_xyz(1:end, 1), trajectory_xyz(1:end, 2));

function showTrajectoryXYZWithoutB(trajectory_xyz)
figure;
scatter3(trajectory_xyz(1:end, 1), trajectory_xyz(1:end, 2), trajectory_xyz(1:end, 3));

function showVxyzAndTime(vxyz, t)
v = sqrt(vxyz(1:end, 1).^2 + vxyz(1:end, 2).^2 + vxyz(1:end, 3).^2);
figure;
plot(t, v);
