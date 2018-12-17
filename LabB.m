%%run LabA.m and get the values first.
%4.5
C = [1 0 0 0];

con = [B A*B A^2*B A^3*B]
obs= [C; C*A;C*A^2;C*A^3]

rankcon = rank(con) %Now: Full rank in con and obs. We've controll
rankobs = rank(obs)  


%%
%4.6
%run LabA.m!
clc

currentPoles = pole(G) %Move integrator and the positive pole %This is TF, might lose some poles/zeros doing like this, CAREFUL
desiredPoles = [currentPoles(2) currentPoles(4) -2 -3]


%adding a feedback K


%syms v K1 K2 K3 K4 b1 b2 b3 b4   %Överger standard sättet för ackerman!
%Asym =sym('a',[4,4])
%Ksym = [K1 K2 K3 K4]
%Bsym = [b1; b2; b3; b4]
%Adet = (Asym-Bsym*Ksym)
%Acl = v*eye(size(A))-(A-B*K)
%Adet = det(Acl) 


%K = [0 0 0 1];

K = acker(A,B,desiredPoles)
disp('from sol')
disp([ -10.0000  -57.4908 -105.0371  -19.5009 ])

%%
%task 4.7
C1 = [8 1 9 3] %try to find some good roots from here

[num den] = ss2tf(A,B,C1,D,1)  %new tf due to new C...
 %[z,p,gain ]  = ss2zp(A, B, C1, D, 1) %find all the zeros/poles!


%find D(s) N(s) and plot Root Locus
s = tf('s');
Ds = den(1)*s^4+den(2)*s^3+den(3)*s^2+den(4)*s+den(5)
Ns =num(1)*s^4+num(2)*s^3+num(3)*s^2+num(4)*s+num(5)

Dsmin = den(1)*s^4-den(2)*s^3+den(3)*s^2-den(4)*s+den(5)
Nsmin = num(1)*s^4-num(2)*s^3+num(3)*s^2-num(4)*s+num(5)
%rlocus(Ns*Nsmin/(Ds*Dsmin))
%axis([-10 10 -5 5])

%more task 4.7
%acker osv
%clc
rho = 5;   %choose from inspection earlier!

% compute the roots of the SRL equation 
afUnsortedAllRoots = rlocus( Ns*Nsmin/(Ds*Dsmin), rho );
[~, aiSortingIndexes] = sort( real(afUnsortedAllRoots) );
afSortedAllRoots = afUnsortedAllRoots(aiSortingIndexes);
afSortedRoots = afSortedAllRoots(1:4)
% compute the gains matrix for the controller K = acker(A, B, afSortedRoots);

% compute the gains matrix for the controller 
K = place(A, B, afSortedRoots)
min(afSortedRoots)
%%
%task 4.8
%full lueberge
C= [1 0 0 0;0 0 1 0];
L = (place(A',C',2*afSortedRoots))'
%%

%task 4.8 reduced Luberger
clc
Cacc = [1 0 0 0];
CNacc = [0 0 1 0];

Vacc = [0 1 0 0;0 0 0 1];
%VNacc = [1 0 0 0;0 1 0 0;0 0 0 1];

invT = [Cacc;CNacc;Vacc]
%invTN = [CNacc;VNacc]
T = inv(invT)
%New basis 1, New sys

AS = invT*A*T
BS = invT*B
BSy = BS(1,1)
BSx = BS(2:4)


CSacc =  Cacc*T
CSNacc = CNacc*T
Cy =CSNacc(1:1)
Cx = CSNacc(2:4)

Ayy =AS(1:1,1:1)
Ayx =AS(1,2:4)
Axy = AS(2:4,1)
Axx = AS(2:4,2:4)

CC = [Ayx;Cx]
L_red=place(Axx',CC',3*afSortedRoots(2:4))' %NO IDEA
L_red_acc = L_red(1:3, 1);
L_red_not_acc = L_red(1:3, 2);


M1 = Axx-L_red_acc*Ayx-L_red_not_acc*Cx
M2 = BSx-L_red_acc*BSy
M3 = Axy-L_red_acc*Ayy-L_red_not_acc*Cy
M4 = L_red_not_acc
M5 = L_red_acc %oklart...
M6 = T(1:4,1)
M7 = T(1:4,2:4)   %only intrested in the xx part???

%L = L_red;

 
%%
%task 4.9
clc
fSamplingPeriod = 0.01;
sys = ss(A,B,C,D)
sys_d = c2d(sys,fSamplingPeriod,'zoh') %sys is sys = ss(A,B,C,D) :) 
[Ad Bd Cd Dd] = ssdata(sys_d)   
%pole(sys_d)
poles_d = exp(afSortedRoots.*fSamplingPeriod)
poles_d_Ld = exp(afSortedRoots.*fSamplingPeriod*3)
Cd = [1 0 0 0;0 0 1 0]
Kd = place(Ad, Bd, poles_d)
Ld = (place(Ad',Cd',poles_d_Ld))'


%see LabB_ObserverOverSimulator_Discrete_Parameters.m 
% :P 



%%
%testa some copy pasta från tidigare!


C1 = [5 1 7 3.5] %try to find some good roots from here
 [num ,den] = ss2zp(Ad,Bd,C1,D,1)
%[num den] = ss2tf(Ad,Bd,C1,D,1)  %new tf due to new C...
%num = real(num)
%[z,p,gain ]  = ss2zp(A, B, C1, D, 1) %find all the zeros/poles!


%find D(s) N(s) and plot Root Locus
s = tf('s');
Ds = den(1)*s^3+den(2)*s^2+den(3)*s+den(4)
Ns =num(1)*s^2+num(2)*s^1+num(3)

Dsmin = -den(1)*s^3+den(2)*s^2-den(3)*s+den(4)
Nsmin =num(1)*s^2-num(2)*s^1+num(3)


%for 5 s'es
%Ds = den(1)*s^4+den(2)*s^3+den(3)*s^2+den(4)*s+den(5)
%Ns =num(1)*s^4+num(2)*s^3+num(3)*s^2+num(4)*s+num(5)
%Dsmin = den(1)*s^4-den(2)*s^3+den(3)*s^2-den(4)*s+den(5)
%Nsmin = num(1)*s^4-num(2)*s^3+num(3)*s^2-num(4)*s+num(5)

rlocus(Ns*Nsmin/(Ds*Dsmin))
%axis([-10 10 -5 5])

%more task 4.7
%acker osv
%clc
rho = 10;   %choose from inspection earlier!

% compute the roots of the SRL equation 
afUnsortedAllRoots = rlocus( Ns*Nsmin/(Ds*Dsmin), rho );
[~, aiSortingIndexes] = sort( real(afUnsortedAllRoots) );
afSortedAllRoots = afUnsortedAllRoots(aiSortingIndexes);
afSortedRoots = afSortedAllRoots(1:4)
% compute the gains matrix for the controller K = acker(A, B, afSortedRoots);

% compute the gains matrix for the controller 
%Kd = place(Ad, Bd, afSortedRoots)
min(afSortedRoots)

%%
%W = C1' * C1;
%Q = rho * W;
%R = 1;
%[Kd, S_lqr, e_lqr] = lqr(Ad, Bd, Q, R) %SNODD KOD.....

%%
C= [1 0 0 0;0 0 1 0];



Caccd = [1 0 0 0];
CNaccd = [0 0 1 0];

Vaccd = [0 1 0 0;0 0 0 1];
%VNacc = [1 0 0 0;0 1 0 0;0 0 0 1];

invTd = [Caccd;CNaccd;Vaccd]
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
%L_redd=place(Axxd',CCd',3*afSortedRoots(2:4))' %NO IDEAx
L_redd = (place(Ad',Cd',3*poles_d))'
L_red_accd = L_redd(1:3, 1);
L_red_not_accd = L_redd(1:3, 2);

            
Md1 = Axxd-L_red_accd*Ayxd-L_red_not_accd*Cxd
Md2 = BSxd-L_red_accd*BSyd
Md3 = Axyd-L_red_accd*Ayyd-L_red_not_accd*Cyd
Md4 = L_red_not_accd
Md5 = L_red_accd %oklart...
Md6 = Td(1:4,1)
Md7 = Td(1:4,2:4)   %only intrested in the x


fGyroConversionFactor  = -1/131;



