g = sym(9.8);
b_f = sym(0);%zero?! 
m_b = sym(0.381); 
l_b = sym(0.112); 
I_b = sym(0.00616);
m_w = sym(0.036);
l_w = sym(0.021);
I_w = sym(0.00000746);
R_m = sym(4.4);
L_m = sym(0);
b_m = sym(0);
K_e = sym(0.444);
K_t = sym(0.470);

%%

%task 3.3.1 SS-form and linerazied
clear
clc
syms g b_f m_b l_b I_b m_w l_w I_w R_m L_m b_m K_e K_t 


gamma = sym('gamma',[2,2])
alpha = sym('alpha',[2,4])
beta = sym('beta',[2,1])
%x = sym('x',[1,4]) %
%syms u
%syms xw xwd xwdd theta thetad thetadd;
%x = [xw xwd theta thetad]'
%xd = [xwdd xwd thetadd thetad]'
%EOM 1
gamma(1,1) =  ( (I_w)/(l_w) + l_w * m_b + l_w * m_w ); 
gamma(1,2) =  m_b * l_b * l_w; 
alpha(1,2) = - ( (K_e * K_t)/(R_m) + b_f ) / (l_w); 
alpha(1,4) =( (K_e * K_t)/(R_m) + b_f );

gamma(2,1) = m_b * l_b;
gamma(2,2) = ( I_b + m_b * l_b^2 ); 
alpha(2,2) =( gamma(2,2) * alpha(1,2) - gamma(1,2) * alpha(2,2)) / delta; 
alpha(2,3) = ( - gamma(1,2) * alpha(2,3)) / delta; 
alpha(2,4) =- ( (K_e * K_t)/(R_m) + b_f ); 







%%
syms a_22 a_23 a_24  a_42 a_43 a_44
Aa=[0 1 0 0; 0  a_22 a_23 a_24;0 0 0 1; 0 a_42 a_43 a_44]
%%
syms g b_f m_b l_b I_b m_w l_w I_w R_m L_m b_m K_e K_t 

gamma_11 = ( (I_w)/(l_w) + l_w * m_b + l_w * m_w ); 
gamma_12 = + m_b * l_b * l_w; 
alpha_11 = 0;
alpha_12 = - ( (K_e * K_t)/(R_m) + b_f ) / (l_w); 
alpha_13 = 0;
alpha_14 = + ( (K_e * K_t)/(R_m) + b_f );
gamma_21 = m_b * l_b;
gamma_22 = ( I_b + m_b * l_b^2 ); 
alpha_22 = ( (K_e * K_t)/(R_m) + b_f ) / (l_w); 
alpha_23 = m_b * l_b * g; 
alpha_24 = - ( (K_e * K_t)/(R_m) + b_f ); 

delta = gamma_11 * gamma_22 - gamma_12 * gamma_21;

beta_11 = + (K_t)/(R_m); beta_12 = + l_w; 
beta_21 = - (K_t)/(R_m); 
beta_22 = + l_b; delta = gamma_11 * gamma_22 - gamma_12 * gamma_21;
a_22 = ( gamma_22 * alpha_12 - gamma_12 * alpha_22) / delta; 
a_23 = ( - gamma_12 * alpha_23) / delta; 
a_24 = ( gamma_22 * alpha_14 - gamma_12 * alpha_24) / delta; 
a_42 = (- gamma_21 * alpha_12 + gamma_11 * alpha_22) / delta; 
a_43 = ( + gamma_11 * alpha_23) / delta;
a_44 = (- gamma_21 * alpha_14 + gamma_11 * alpha_24) / delta; 
b_21 = ( gamma_22 * beta_11 - gamma_12 * beta_21) / delta; 
b_22 = ( gamma_22 * beta_12 - gamma_12 * beta_22) / delta; 
b_41 = (- gamma_21 * beta_11 + gamma_11 * beta_21) / delta;
b_42 = (- gamma_21 * beta_12 + gamma_11 * beta_22) / delta; 

%
A = [ 0 1 0 0 ; 0 a_22 a_23 a_24 ; 0 0 0 1 ; 0 a_42 a_43 a_44 ];
B = [ 0 ; b_21 ; 0 ; b_41 ]; 
%C = [ 0 0 1 0 ]; 

% 
% B with an additional input for the poking disturbances 
Bd = [ 0 0 ; b_21 b_22 ; 0 0 ; b_41 b_42 ]; 





%R1 =( inv(gamma)*alpha)*x
%R2 = inv(gamma)*beta*u
%xd = R1+R2