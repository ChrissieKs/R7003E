%Run LabA Firstly.
%% Task 5.1.1 Last thing to do
freq = 200;                 %FIND BEST FREQ HERE :) 
fSamplingPeriod = 1/freq;

sys = ss(A,B,C,D)
sys_d = c2d(sys,fSamplingPeriod,'zoh') 
[Ad Bd Cd Dd] = ssdata(sys_d)  
Cd = [1 0 0 0;0 0 1 0];

%% Task 5.2.1 Discrete LQR
%For LQR(A,B,Q,R).. Q = C'*C*rho.... R=1 ? 
rho = 1000;
Q = rho*Cd'*Cd;
R =1;
Kd = dlqr(Ad,Bd,Q,R)
%For rho = 1, smallest freq = 5


%% Task 5.3.1
%pz' = e^(x)*pz  moving poles.
%2
freq =200;                  
fSamplingPeriod = 1/freq+0.00;

sys = ss(A,B,C,D)
sys_d = c2d(sys,fSamplingPeriod,'zoh') 
[Ad Bd Cd Dd] = ssdata(sys_d)  
Cd = [1 0 0 0;0 0 1 0];

Acroots = eig(Ad-Bd*Kd) %roots for closed???
%3 We need to find roots
%Using 2nd order? Is it even possible here? 
%LQR/Root locus? - Let's try.

%Kd = place(Ad, Bd, Acroots)
Ld= (place(Ad',Cd',Acroots.^3))'   

%
%
%
%
% Reduced
%using eig as (Ad-Bd*Kd) one probably should cahnge Ad = Ad-Bd*Kd
C= [1 0 0 0;0 0 1 0];



Caccd = [1 0 0 0];
CNaccd = [0 0 1 0]; 

Vaccd = [0 1 0 0;0 0 1 0;0 0 0 1]; 
%VNacc = [1 0 0 0;0 1 0 0;0 0 0 1];

invTd = [Caccd;Vaccd]
%invTN = [CNacc;VNacc]
Td = inv(invTd)
%New basis 1, New sys

ASd = invTd*Ad*Td
BSd = invTd*Bd
BSyd = BSd(1,1)
BSxd = BSd(2:4)


CSaccd =  Caccd*Td
CSNaccd = CNaccd*Td
Cyd =CSNaccd(1:1)
Cxd = CSNaccd(2:4)

Ayyd =ASd(1:1,1:1)
Ayxd =ASd(1,2:4)
Axyd = ASd(2:4,1)
Axxd = ASd(2:4,2:4)

CCd = [Ayxd;Cxd]
afpoles_d_Ld = sort(real(Acroots));
L_redd=place(Axxd',CCd',afpoles_d_Ld(2:4))' %Probably issue here!
%L_redd = Ld%(place(Ad',Cd',3*poles_d))'
L_red_accd = L_redd(1:3, 1);
L_red_not_accd = L_redd(1:3, 2);

            
Md1 = Axxd-L_red_accd*Ayxd-L_red_not_accd*Cxd
Md2 = BSxd-L_red_accd*BSyd
Md3 = Axyd-L_red_accd*Ayyd-L_red_not_accd*Cyd
Md4 = L_red_not_accd
Md5 = L_red_accd %oklart...
Md6 = Td(1:4,1)
Md7 = Td(1:4,2:4)   %only intrested in the x

%% 5.5 External Reference
Cd_x_w = [1 0 0 0]; %oklart....känns som Cd behövs ändras...

MdN = acker((Ad-Bd*Kd-Ld*Cd)',Kd',Acroots)'

%Cd = Cd_x_w
Cd1 = Cd_x_w
%ao = 
%ac = 
test1 = Cd1*Bd*Kd*Bd
test2 = Cd1*(Ad-Ld*Cd)*Bd
test3 =ao*ac;   
testINV = inv(Cd1*Bd*Kd)

M = testINV*(test1-test2+test3) %Derivation on paper.
%M/N = acker....
%N = M/acker...=M/(MdN)
N = M'./MdN
Nud = N
Nxd =0;


%% 5.5 External 2nd Try
MdN = acker((Ad-Bd*Kd-Ld*Cd)',Kd',Acroots)'
Cd_x_w = [1 0 0 0]; %oklart....känns som Cd behövs ändras...

Nud = -inv(Cd_x_w*inv(Ad-Bd*Kd)*Bd)
Nxd = 0;









