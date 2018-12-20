%kör labB först obviously.
%% 4.5
clc
A = sym('a', [4 4])
B = sym('b',[4 1]) 
%jaaa...
syms A B
[B A*B A^2*B A^3*B]
%% 4.6
afFigurePosition = [1 1 10 6];

figure(2)
plot(theta_b_4_6.time, theta_b_4_6.signals.values * 180 / pi);
title('\theta_b');
xlabel('time');
ylabel('degrees')
legend({'\theta_b','lin:?_b'},'Location','southeast')
print('-depsc2', '-r300', 'C:\Users\Oscar\Desktop\test\MatlabAndSimulink\Bilder\Task_4_6_theta_b.eps');


figure(3)
plot(u.time, u.signals.values); 
title('v_m');
xlabel('time');
ylabel('Volt') 
set(gcf, 'Units', 'centimeters'); 
set(gcf,'Position',afFigurePosition); 
set(gcf, 'PaperPositionMode', 'auto');
legend({'v_m','lin:v_m'},'Location','northwest')
print('-depsc2', '-r300', 'C:\Users\Oscar\Desktop\test\MatlabAndSimulink\Bilder\Task_4_6_v.eps');
%% 4.7
figure(2)
plot(theta_b_4_6.time, theta_b_4_6.signals.values * 180 / pi);
title('\theta_b');
xlabel('time');
ylabel('degrees')
legend({'\theta_b','lin:?_b'},'Location','southeast')
print('-depsc2', '-r300', 'C:\Users\Oscar\Desktop\test\MatlabAndSimulink\Bilder\Task_4_7_theta_b.eps');


figure(3)
plot(u.time, u.signals.values); 
title('v_m');
xlabel('time');
ylabel('Volt') 
set(gcf, 'Units', 'centimeters'); 
set(gcf,'Position',afFigurePosition); 
set(gcf, 'PaperPositionMode', 'auto');
legend({'v_m','lin:v_m'},'Location','northeast')
print('-depsc2', '-r300', 'C:\Users\Oscar\Desktop\test\MatlabAndSimulink\Bilder\Task_4_7_v.eps');

%% 4.8

figure(2)
plot(theta_b.time, theta_b.signals.values * 180 / pi);
title('\theta_b');
xlabel('time');
ylabel('degrees')
legend({'\theta_b of PID/"Simulation"','\theta_b of Full Order','\theta_b of Reduced Order'},'Location','northeast')
print('-depsc2', '-r300', 'C:\Users\Oscar\Desktop\test\MatlabAndSimulink\Bilder\Task_4_8_theta_b.eps');

figure(20)
plot(x_w.time, x_w.signals.values);
title('\x_w');
xlabel('time');
ylabel('Distance [m]')
legend({'x_w','full x_w','reduced x_w'},'Location','southeast')
print('-depsc2', '-r300', 'C:\Users\Oscar\Desktop\test\MatlabAndSimulink\Bilder\Task_4_8_x_w.eps');

%% 4.8 Error!
theta_b_error_full = abs(theta_b.signals.values(1:end,2) - theta_b.signals.values(1:end,1))* 180 / pi

theta_b_error_red = abs(theta_b.signals.values(1:end,3) - theta_b.signals.values(1:end,1))* 180 / pi

max_error_full = max(theta_b_error_full)

max_error_red = max(theta_b_error_red)



%% 4.9

figure(491)
stairs(theta_b.time, theta_b.signals.values(1:end,1) * 180 / pi);
title('\theta_b Discrete with full Luenberger ');
xlabel('time');
ylabel('degrees')
legend({'\theta_b'},'Location','northeast')
print('-depsc2', '-r300', 'C:\Users\Oscar\Desktop\test\MatlabAndSimulink\Bilder\Task_4_9_theta_b.eps');

figure(492)
plot(x_w.time, x_w.signals.values(1:end,1));
title('x_w Discrete with full Luenberger ');
xlabel('time');
ylabel('Distance [m]')
legend({'x_w'},'Location','northeast')
print('-depsc2', '-r300', 'C:\Users\Oscar\Desktop\test\MatlabAndSimulink\Bilder\Task_4_9_x_w.eps');



figure(493)
plot(u.time, u.signals.values); 
title('v_m Discrete with full Luenberger ');
xlabel('time');
ylabel('Volt') 
set(gcf, 'Units', 'centimeters'); 
set(gcf,'Position',afFigurePosition); 
set(gcf, 'PaperPositionMode', 'auto');  
legend({'v_m','lin:v_m'},'Location','northeast')
print('-depsc2', '-r300', 'C:\Users\Oscar\Desktop\test\MatlabAndSimulink\Bilder\Task_4_9_v.eps');
















