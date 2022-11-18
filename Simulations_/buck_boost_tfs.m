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
Cb = Vo*T*D/2/Rb/delta_vo
Lb = Vs*D*T/2/delta_iL
% Results at D duty cyclce approximately confirm the ripple voltage. (constant D not PI)
% Performance degrades above 0.5 duty cycle 

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

%% PI design
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
