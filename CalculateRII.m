function CalculateRII

syms theta r phi phi_1 c1 c2 c3 a w t 

cross([1, 0, 0], cross([1, 0, 0], [c1, c2 c3]))

a = [sin(theta) * cos(phi) sin(theta) * sin(phi) cos(theta)];

b = [cos(phi_1) sin(phi_1), 0];

func = cos(phi) * exp(i * phi)
int(func)

getVectorAngle(a, b)

int(cos(w * t - phi), t)

%f = int(exp(-j * phi_1) / (c2 - c1 * cos(phi - phi_1) + c3 * cos(phi - phi_1) * cos(phi - phi_1)));
%f0 = subs(f, phi_1, 0);
%f2pi = subs(f, phi_1, 2 * pi);
%pretty(f)

syms a
v_dot = [-a * cos(phi_1) -a * sin(phi_1) 0];
en = [sin(theta) * cos(phi) sin(theta) * sin(phi) cos(theta)];
v = cross(en, cross(en, v_dot))

getVectorAngle(v_dot, en);
sin(getVectorAngle(en, cross(en, v_dot))) * sin(getVectorAngle(v_dot, en));


x = 0:0.1:20
plot(x, besselj(1,sin(x)))
figure
plot(x, besselj(-1,x))
f = int(exp(1i * (phi)) * exp(1i * c1 * cos(phi)), phi,0, 2 * pi)

1i * exp(-1i * phi) 


x = 0 : 0.1 : 20
y = 0 : 0.1 : 20
f = linspace(0, 2 * pi, 200)
an = i * exp(i * f) - 1 / i * exp(-i * f) 
uv = exp(i * f) - exp(-i * f) 
plot(angle(uv))

function valueabs = getValueAbs(value)
valueabs = sqrt(sum(value .* value));

function angle = getVectorAngle(a, b)
c=dot(a,b);  % 求内积
 
d=dot(a,a);  % 求a的长度
e=sqrt(d);
 
f=dot(b,b) ; % 求b的长度
g=sqrt(f);
 
h=c/(e*g);
 
z = acos(h);
angle = z;

