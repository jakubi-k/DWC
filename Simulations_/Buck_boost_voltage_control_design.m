%% System parameters
Vs = 415;
Vb = 400;
s_fre = 85e3;
T = 1/s_fre;
Rb = 1;

%% Target performance
Vo = 412;
delta_vo = 0.5; %Ripple voltage
zeta = 1/sqrt(2);

%% Inductor/Capacitor values
D = Vo/(Vs+Vo);
Cb = (Vo-Vb)*T/2/Rb/delta_vo;
Lb =(2*Rb*D*zeta)^2*Cb;
% Results at D duty cyclce approximately confirm the ripple voltage. (constant D not PI)
% Performance degrades above 0.5 duty cycle 

s = tf('s');


%%
kp = 4;
ki = 10;
Gc =  kp + ki/s;
T = 1/(1/Gc/Gp+1);
step(T)
%%

% tr = 25e-1;
% wn = 4.25/zeta/tr;
% zeta = 1/sqrt(2);
% a = Lb*Cb;
% b = Lb/Rb;
% c = D^2;
% 
% alpha = b/a-2*zeta*wn;
% ki = a*alpha*wn^2
% kp = (2*alpha*zeta*wn+wn^2)*a-c
% 
% Gc =  kp + ki/s;
% T = 1/(1/Gc/Gp+1);
% figure(1);step(T)
% figure(2);hold on;pzmap(T)

% theory apraoch above needs checking - currently not ideal
ki = 1e3;
kp = 1;
T = 1/(1/Gc/Gp+1);
figure(1);step(T)
figure(2);pzmap(T)