%% System parameters
Vs = 415;
Vb = 0;
s_fre = 85e3;
T = 1/s_fre;
Rb = 1;

%% Target performance
D = 0.4;
Vo = Vs*D/(1-D);
delta_vo = 0.5; %Ripple voltage
delta_iL = 0.5;
zeta = 1/sqrt(2);

%% Inductor/Capacitor values
Cb = Vo*T*D/2/Rb/delta_vo;
Lb = Vs*D*T/2/delta_iL;
% Results at D duty cyclce approximately confirm the ripple voltage. (constant D not PI)
% Performance degrades above 0.5 duty cycle 

%%
s = tf('s');
Gp = -Vs*(s*Lb*D/Rb/(1-D)^2-1)/(s^2*Lb*Cb+s*Lb/Rb+(1-D)^2)/Rb;
kp = 0.0005;
ki = 0.1;
Gc = (s*kp+ki)/s;
T = 1/(1+1/Gp/Gc);

% num = Vs*(ki + kp*s)*(Rb + D^2*Rb - 2*D*Rb - D*Lb*s);
% den = D^4*Rb^2*s - 4*D^3*Rb^2*s + Cb*Lb*D^2*Rb^2*s^3 + 6*D^2*Rb^2*s + Lb*D^2*Rb*s^2 + Vs*kp*D^2*Rb*s + Vs*ki*D^2*Rb - 2*Cb*Lb*D*Rb^2*s^3 - 4*D*Rb^2*s - 2*Lb*D*Rb*s^2 - 2*Vs*kp*D*Rb*s - 2*Vs*ki*D*Rb - Lb*Vs*kp*D*s^2 - Lb*Vs*ki*D*s + Cb*Lb*Rb^2*s^3 + Rb^2*s + Lb*Rb*s^2 + Vs*kp*Rb*s + Vs*ki*Rb;

a = Lb*Cb*Rb*(1-D)^2;
b = Lb*Rb*(1-D)^2-Lb*Vs*kp*D;
c = Vs*kp*Rb*(1-D)^2*s + (1-D)^4*Rb^2*s - Lb*Vs*ki*D*s ;
d = Vs*ki*Rb*(1-D)^2;

T2 = 1/(a*s^3+b*s^2+c*s+d);


% figure; pzmap(Gp);
figure(1)
pzmap(T);
figure; step(T)
xlim([0, 0.1])
ylim([-0.2, 1.2])

% hold on;
% pzmap(T2);
%% TF
s = tf('s');
Vo = Vs*D/(1-D);
IL=Vo/Rb/(1-D);
Gp = -Vs*(s*Lb*D/Rb/(1-D)^2-1)/(s^2*Lb*Cb+s*Lb/Rb+(1-D)^2);

figure
subplot(2,1,1)
plot(io.Time,io.Data)
ylabel("Amplitude")
xlim([0.04 0.06])
grid on
% ylim([0.75 1.2])
title('Simulation results')
hold on
subplot(2,1,2)
step(Gp)
grid on
xlim([0 0.02])
% ylim([0 1.2])
title('transfer function:step response')

%% PI design - analysis - incorrect
syms vs s D Lb Rb Cb D_ kp ki Vs
Gp = -Vs*(s*Lb*D/Rb/(1-D)^2-1)/(s^2*Lb*Cb+s*Lb/Rb+(1-D)^2)/Rb;
Gc = (s*kp+ki)/s;
T = 1/(1+1/Gp/Gc);
ex=expand(T);
T = simplify(ex);
[num, den] = numden(T);

syms zeta wn alpha

a = Lb*Cb*Rb*(1-D)^2;
b = Lb*Rb*(1-D)^2-Lb*Vs*kp*D;
c = Vs*kp*Rb*(1-D)^2*s + (1-D)^4*Rb^2*s - Lb*Vs*ki*D*s ;
d = Vs*ki*Rb*(1-D)^2;

eq1 = 2*zeta*wn+alpha == b/a;
eq2 = wn^2+2*alpha*wn+alpha*wn^2 == c/a;
eq3 = alpha*wn^2 == d/a;
zwa = solve([eq1, eq2, eq3], [zeta, wn, alpha])
%% Kp Ki
