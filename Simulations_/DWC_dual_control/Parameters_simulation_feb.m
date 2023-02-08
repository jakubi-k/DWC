speed=20;%100*60^2/1.609e3;%check mph or ms
a1 = 0.3;
b1 = 0.3;
a2 = 0.3;
b2 = 0.3;

h = 0.3; % distance between coils (air gap)


%MOSFET
Rds=80e-3;
%PWM
m = 0.9; %modulation index
switching_ratio = 27;
%parameters
% f = 50;
phase_a = 0;
%%% Iteration
%% Preliminary parameters
V_load = 300/m; %tesla model S/X batteries 300V(75kWh)-400V(100kWh)
P_load = 3.7e3; 
% P_load = 10e3; 
R_load = (V_load^2)/P_load;
Vp = 300;
Kw = 1;
f_opt = 20e3;
d1_max = 4*10; %[A/mm2]
d2_max = 4*10; %[A/mm2]
N1_max = 25*100;
N2_max = 25*100;
mu0 = 4*pi*1e-7;
rho_cu = 17.2e-9; %resistivity copper annealed ______ DON'T FORGET TO CHANGE IT!!!!!!!!!
P2 = 0;
V2 = 0;
%% transmitter coil
N1 = N1_max;
S1 = 1e-6;
% a1 = 0.1;
% b1 = 0.1;

%% receiver coil
N2 = N2_max;
S2 = 1e-6;
% a2 = 0.1;
% b2 = 0.1;

%% MUTUAL INDUCTANCE
c = 0.00; % x distance from z axis
e = 0.00; % y distance from z axis

%change of variable
d = a1 -a2 - c;
m = a2 + c;
q = b2 + e;
g = a1 - c;
p = b1 - e;
t = b1 - b2 - e;

while 1
    %Transmitter coil parameters
    R1 = rho_cu*N1*2*(a1+b1)/S1;
    r1 = sqrt(N1*S1/pi); %equivalent radius of the winding
    L1 = (mu0/pi)*(N1^2)*(a1*log((2*a1*b1)/(r1*(a1+sqrt(a1^2 + b1^2)))) + b1*log((2*a1*b1)/(r1*(b1+sqrt(a1^2 + b1^2)))) - 2*(a1 +b1 - sqrt(a1^2 + b1^2)) + 0.25*(a1 + b1));
    %Receiver coil parameters
    R2 = rho_cu*N2*2*(a2+b2)/S2;
    r2 = sqrt(N2*S2/pi); %equivalent radius of the winding
    L2 = (mu0/pi)*(N2^2)*(a2*log((2*a2*b2)/(r2*(a2+sqrt(a2^2 + b2^2)))) + b2*log((2*a2*b2)/(r2*(b2+sqrt(a2^2 + b2^2)))) - 2*(a2 +b2 - sqrt(a2^2 + b2^2)) + 0.25*(a2 + b2));
    %Mutual inductance calculation
    M1 = ((d*log((d+sqrt(h^2 + (-t)^2 + d^2))/(d + sqrt(h^2 + d^2 + q^2)))) + (g*log((g+sqrt(h^2 + (q)^2 + g^2))/(g + sqrt(h^2 + g^2 + (-t)^2)))) + (c*log(((-c)+sqrt(h^2 + (q)^2 + c^2))/((-c) + sqrt(h^2 + c^2 + (-t)^2)))) + (m*log(((-m)+sqrt(h^2 + (-t)^2 + m^2))/((-m) + sqrt(h^2 + m^2 + q^2)))) + sqrt(h^2 + q^2 + d^2) - sqrt(h^2 + q^2 + g^2) - sqrt(h^2 + q^2 + m^2) + sqrt(h^2 + q^2 + c^2) + sqrt(h^2 + (-t)^2 + g^2) - sqrt(h^2 + (-t)^2 + d^2) + sqrt(h^2 + (-t)^2 + m^2) - sqrt(h^2 + (-t)^2 + c^2));
    M2 = ((d*log((d+sqrt(h^2 + (-p)^2 + d^2))/(d + sqrt(h^2 + d^2 + e^2)))) + (g*log((g+sqrt(h^2 + (e)^2 + g^2))/(g + sqrt(h^2 + g^2 + (-p)^2)))) + (c*log(((-c)+sqrt(h^2 + (e)^2 + c^2))/((-c) + sqrt(h^2 + c^2 + (-p)^2)))) + (m*log(((-m)+sqrt(h^2 + (-p)^2 + m^2))/((-m) + sqrt(h^2 + m^2 + e^2)))) + sqrt(h^2 + e^2 + d^2) - sqrt(h^2 + e^2 + g^2) - sqrt(h^2 + e^2 + m^2) + sqrt(h^2 + e^2 + c^2) + sqrt(h^2 + (-p)^2 + g^2) - sqrt(h^2 + (-p)^2 + d^2) + sqrt(h^2 + (-p)^2 + m^2) - sqrt(h^2 + (-p)^2 + c^2));
    M3 = ((t*log((t+sqrt(h^2 + (-g)^2 + t^2))/(t + sqrt(h^2 + t^2 + c^2)))) + (p*log((p+sqrt(h^2 + (p)^2 + c^2))/(p + sqrt(h^2 + (-g)^2 + p^2)))) + (e*log(((-e)+sqrt(h^2 + (e)^2 + c^2))/((-e) + sqrt(h^2 + e^2 + (-g)^2)))) + (q*log(((-q)+sqrt(h^2 + (-g)^2 + q^2))/((-q) + sqrt(h^2 + c^2 + q^2)))) + sqrt(h^2 + c^2 + t^2) - sqrt(h^2 + c^2 + p^2) - sqrt(h^2 + c^2 + q^2) + sqrt(h^2 + e^2 + c^2) + sqrt(h^2 + (-g)^2 + p^2) - sqrt(h^2 + (-g)^2 + t^2) + sqrt(h^2 + (-g)^2 + q^2) - sqrt(h^2 + (-g)^2 + e^2));
    M4 = ((t*log((t+sqrt(h^2 + (-d)^2 + t^2))/(t + sqrt(h^2 + t^2 + m^2)))) + (p*log((p+sqrt(h^2 + (m)^2 + p^2))/(p + sqrt(h^2 + (-d)^2 + p^2)))) + (e*log(((-e)+sqrt(h^2 + (e)^2 + m^2))/((-e) + sqrt(h^2 + e^2 + (-d)^2)))) + (q*log(((-q)+sqrt(h^2 + (-d)^2 + q^2))/((-q) + sqrt(h^2 + m^2 + q^2)))) + sqrt(h^2 + m^2 + t^2) - sqrt(h^2 + m^2 + p^2) - sqrt(h^2 + m^2 + q^2) + sqrt(h^2 + e^2 + m^2) + sqrt(h^2 + (-d)^2 + p^2) - sqrt(h^2 + (-d)^2 + t^2) + sqrt(h^2 + (-d)^2 + q^2) - sqrt(h^2 + (-d)^2 + e^2)); 
    M = (mu0/(4*pi))*N1*N2*(-M1 + M2 + M3 - M4);
    %Calculations
    w_opt = Kw*sqrt(R1*(R2+R_load))/M;
    f_opt = w_opt/(2*pi);
    C1 = 1/((w_opt^2)*L1);
    C2 = 1/((w_opt^2)*L2);
    I1 = Vp/(R1 + w_opt*L1*1i - (1/(w_opt*C1))*1i + ((w_opt * M)^2)/(R2 + R_load + w_opt*L2*1i - (1/(w_opt*C2))*1i));
    I1_abs = abs(I1);
    I2 = ((-w_opt*M*1i)/(R2 + R_load + w_opt*L2*1i - (1/(w_opt*C2))*1i))*I1;
    I2_abs = abs(I2);
    V2 = I2_abs*R_load;
    V_C1 = I1*(-1i/(w_opt*C1));
    V1_C1_abs = abs(V_C1);
    V_C2 = I2*(-1i/(w_opt*C2));
    V2_C2_abs = abs(V_C2);
    d1 = I1_abs/(S1*1e6);
    d2 = I2_abs/(S2*1e6);
    Qp = (L1*R_load)/(w_opt*M^2);
    Qs = (w_opt*L2)/R_load;
    P2 = R_load*(I2_abs^2);
    n = ((R_load*(I2_abs^2))/(R1*(I1_abs^2) + R2*(I2_abs^2) + R_load*(I2_abs^2)));
    
    P2
    if (P2 <= P_load)
        d1
        if (d1 < d1_max)
            d2
            if (d2 < d2_max)
               if (((P_load -50)<= P2)||(P2 <= P_load+50)) && (((V_load -5)<= V2)||(V2 <= V_load+5)) && (81.39e3 <= f_opt)&&(f_opt <= 90e3) && (Qp > Qs) && (d1 < d1_max) && (d2 <= d2_max)
                        break;
               else
                   N1
                   if (N1 > 0)
                        N1 = N1 - 1;
                        N2
                        if (N2 > 0)
                            N2 = N2 - 1;
                        else
                            N2 = N2_max;
                        end
                   else
                       N1 = N1_max;   
                   end
               end
            else
                S2 = S2 + 2.5e-7;
            end
        else
            S1 = S1 + 2.5e-7;
        end
    else
       Kw = Kw + 0.01; 
    end
end 
        
%%Switchig frequency
Fc = switching_ratio*f_opt; %Carrier frequency
N1, N2, Fc