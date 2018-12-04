%%
%task 3.3.1 SS-form and linerazied
clear
clc

g = 9.8;
b_f = 0; 
m_b = 0.381; 
l_b = 0.112;   
I_b = 0.00616;
m_w = 0.036;
l_w = 0.021;
I_w = 0.00000746;
R_m = 4.4;
L_m = 0;
b_m = 0;
K_e = 0.444;
K_t = 0.470;

%a = sym('a', [2 4]);
%gamma = sym('gamma',[2,2]);
%beta = sym('beta',[2,1]);
%copy pastad kod


gamma_11 = ( (I_w)/(l_w) + l_w * m_b + l_w * m_w ); gamma_12 = + m_b * l_b * l_w; alpha_12 = - ( (K_e * K_t)/(R_m) + b_f ) / (l_w); alpha_14 = + ( (K_e * K_t)/(R_m) + b_f ); gamma_21 = m_b * l_b; gamma_22 = ( I_b + m_b * l_b^2 ); alpha_22 = ( (K_e * K_t)/(R_m) + b_f ) / (l_w); alpha_23 = m_b * l_b * g; alpha_24 = - ( (K_e * K_t)/(R_m) + b_f ); beta_11 = + (K_t)/(R_m); beta_12 = + l_w; beta_21 = - (K_t)/(R_m); beta_22 = + l_b; delta = gamma_11 * gamma_22 - gamma_12 * gamma_21; a_22 = ( gamma_22 * alpha_12 - gamma_12 * alpha_22) / delta; a_23 = ( - gamma_12 * alpha_23) / delta; a_24 = ( gamma_22 * alpha_14 - gamma_12 * alpha_24) / delta; a_42 = (- gamma_21 * alpha_12 + gamma_11 * alpha_22) / delta; a_43 = ( + gamma_11 * alpha_23) / delta;
a_44 = (- gamma_21 * alpha_14 + gamma_11 * alpha_24) / delta; 
b_21 = ( gamma_22 * beta_11 - gamma_12 * beta_21) / delta; 
b_22 = ( gamma_22 * beta_12 - gamma_12 * beta_22) / delta; 
b_41 = (- gamma_21 * beta_11 + gamma_11 * beta_21) / delta;
b_42 = (- gamma_21 * beta_12 + gamma_11 * beta_22) / delta; 
%
A = [ 0 1 0 0 ; 0 a_22 a_23 a_24 ; 0 0 0 1 ; 0 a_42 a_43 a_44 ];
B = [ 0 ; b_21 ; 0 ; b_41 ]; 
C = [ 0 0 1 0 ]; 
D=0;
% 
% B with an additional input for the poking disturbances 
Bd = [ 0 0 ; b_21 b_22 ; 0 0 ; b_41 b_42 ];  %try to prove this in 3.6.1





%%
%3.4.1
s = tf('s')
sys  = ss(A,B,C,D) %state space
poles = pole(sys) % 4 poles, 1-2(?) unstable
zerosG = zero(sys) %2 zeros at z = 0;
G = s^2/(s*(s-poles(1))*(s-poles(2))*(s-poles(3)))
[a,b,gain] = ss2zp(A,B,C,D,1) %with gain from here
G = G*-90 %-90 = gain...
%[num den ]=ss2tf(A,B,C,D,1)





%ss2tf(A,B,C,D,1)

%%
%3.5.1
p1 = poles(1)
p2 = poles(2)
p3 = poles(3) %den instabila polen :( 

%choose new poles here!
p2_new = p2;
p3_new = p2;    

%can change this value for different poles! will probably be unstable when far away from orignial pole?! 
%enligt eq (58)
%To change 1 POLE, use this!
%p = (p1*p3_new+p2*p3_new-p1*p3-p2*p3)/gain 
%i = (p1*p2*p3-p1*p2*p3_new)/gain
%d = (p3-p3_new)/gain

%to change 2 POLES use this!

p =(p1*p2_new-p1*p2+p1*p3_new-p1*p3+p2_new*p3_new-p2*p3)/gain 
i = (p1*p2*p3-p1*p2_new*p3_new)/gain
d = (p3-p3_new+p2-p2_new)/gain
d_temp = d;

%new TF 
%pidC = 
%G_new_poles = 
%G_pid = ga




%3.5.1
%syms p1 p2 p3 kp ki kd gain s
%s = tf('s')
%Gtemp =  gain*s^2/((s-p1)*(s-p2)*(s-p3))
%PIDtemp = kp+ki/s+kd*s
%feedback(Gtemp,PIDtemp)
%feedbacktemp =  Gtemp*PIDtemp/(1+PIDtemp*Gtemp)



%3.5.1 
%Just playing with some code
%p = -46.6;
%i = -260.4;
%d = -0.1;
Gclp = p*G/(1+p*G) % feedback with P-controller
P_poles = pole(Gclp)
P_zeros=zero(Gclp)
rlocus(Gclp)
%How to do a root locus as in instructions?!

%PI-controller

pidi = pid(p,i,0,0)
Gclpi = feedback(G,pidi) %feedback with PI-controller
PI_poles = pole(Gclpi)
PI_zeros=zero(Gclpi)

pidi = pid(p,i,0,0)
Gclpi = feedback(G,pidi)
PI_poles = pole(Gclpi)
PI_zeros=zero(Gclpi)

%pid-controller
pidid = pid(p,i,d,0)
Gclpid = feedback(G,pidid) %feedback with PI-controller
PID_poles = pole(Gclpid)
PID_zeros=zero(Gclpid)
impulse(Gclpid)
title('Impulse Response for PID Closed loop system')
%rlocus(G*pidid)




%%
%3.4.1 COPY PASTAD NOT IN USE
%s(s-a22)(s^2-s*a44-a42) roots for TF
%clc
%poles = [0 a_22 a_44/2+sqrt(a_44^2/4+a_42) a_44/2-sqrt(a_44^2/4+a_42)]
%zeros = [ 0 a_22 B(4)]
%ppoles doesn't seem to match with solutions..

%tLinearizedBot = ss(A, B, C, 0)
%[tLinearizedBotNumerator, tLinearizedBotDenominator] = ss2tf(A, B, C, 0, 1) 
%[afLinearizedBotZeros, afLinearizedBotPoles, fLinearizedBotGain] = ss2zp(A, B, C, 0, 1)


%%
%plotting 3.4.1
figure(1) 
pzplot(sys); 
grid;
print('-depsc2', '-r300', 'poles_zeros_map_of_linearized_balancing_robot.eps');


%%
%3.6.1
%ta fram Bd.









%%
%task 3.7.1
%syms z
%f = 6*4;
%T = 1/f;
%pididZ = p+i*T*z/(z-1)^2+d*(z-1)/z  %högst oklart
%sysDd = c2d(pidid/s,T,'zoh')            %högst oklart...



C = eye(4)
 D= zeros(4,2)
%%
%plotting stuff
name = ('task_3_7')
open_system('C:\Users\OscarBäckman\Documents\MATLAB\Reglerteknik\task_3_7'); 
saveas(get_param(name,'Handle'),'task_3_7.pdf');
sim(name);
close_system(name);

afFigurePosition = [1 1 10 6];
figure(1) 
plot(x_w.time, x_w.signals.values); 
title('x_w'); xlabel('time');
ylabel('meters')
set(gcf, 'units', 'centimeters'); 
set(gcf,'Position',afFigurePosition); 
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc2', '-r300', 'Task_3_7_x_w.eps');

figure(2)
plot(theta_b.time, theta_b.signals.values * 180 / pi);
title('\theta_b');
xlabel('time');
ylabel('degrees')

set(gcf, 'Units', 'centimeters');
set(gcf,'Position',afFigurePosition);
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc2', '-r300', 'Task_3_7_theta_b.eps');

figure(3) 
plot(d.time, d.signals.values);
title('d'); 
xlabel('time');
ylabel('Newton')
set(gcf, 'Units', 'centimeters');
set(gcf,'Position',afFigurePosition);
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc2', '-r300', 'Task_3_7_d.eps');

figure(4)
plot(v_m.time, v_m.signals.values); 
title('v_m');
xlabel('time');
ylabel('Volt') 
set(gcf, 'Units', 'centimeters'); 
set(gcf,'Position',afFigurePosition); 
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc2', '-r300', 'Task_3_7_v_m.eps'); 


%%
%task 3.8.1
bode(G) % this? %wmax = 5.9
f3db = 2.2282;
f_sample = f3db*25;
T = 1/f_sample
%bode(Gclpid)
%f = 10; %sampling rate
%T = 1/f;
%sysDd = c2d(Gclpid,T,'zoh')


%bode(Gclpid*s) %ehm?

%%
%spara från simunlink saker för task 3.8.1


name = ('task_3_8_b')
open_system('C:\Users\OscarBäckman\Documents\MATLAB\Reglerteknik\task_3_8_b'); 
saveas(get_param(name,'Handle'),'task_3_8_b.pdf');
sim(name);
%close_system(name);

afFigurePosition = [1 1 10 6];

figure(2)
plot(theta_b.time, theta_b.signals.values * 180 / pi);
title('\theta_b');
xlabel('time');
ylabel('degrees')
legend({'?_b','lin:?_b'},'Location','northwest')
print('-depsc2', '-r300', 'Task_3_8_continous_theta_b.eps'); 


figure(4)
plot(v_m.time, v_m.signals.values); 
title('v_m');
xlabel('time');
ylabel('Volt') 
set(gcf, 'Units', 'centimeters'); 
set(gcf,'Position',afFigurePosition); 
set(gcf, 'PaperPositionMode', 'auto');
legend({'v_m','lin:v_m'},'Location','northwest')
print('-depsc2', '-r300', 'Task_3_8_continous_v_m.eps'); 
%%
name = ('task_3_8_discreate');

open_system('C:\Users\OscarBäckman\Documents\MATLAB\Reglerteknik\task_3_8_b'); 
saveas(get_param(name,'Handle'),'task_3_8_dicreate.pdf');
sim(name);
afFigurePosition = [1 1 10 6];

figure(2)
plot(theta_b.time, theta_b.signals.values * 180 / pi);
title('\theta_b');
xlabel('time');
ylabel('degrees')
legend({'?_b','lin:?_b'},'Location','northwest')
print('-depsc2', '-r300', 'Task_3_8_discreate_theta_b.eps'); 


figure(4)
plot(v_m.time, v_m.signals.values); 
title('v_m');
xlabel('time');
ylabel('Volt') 
set(gcf, 'Units', 'centimeters'); 
set(gcf,'Position',afFigurePosition); 
set(gcf, 'PaperPositionMode', 'auto');
legend({'v_m','lin:v_m'},'Location','northwest')
print('-depsc2', '-r300', 'Task_3_8_discreate_v_m.eps'); 


