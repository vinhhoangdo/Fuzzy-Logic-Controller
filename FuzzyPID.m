clear all, close all, clc
%% New FIS GUI
%Create new fis
sys = newfis('SugFuzzyPID','FISType','sugeno')
%% Inputs
% Define input variable error "E" to fis 
%with range [-2 4]

sys = addInput(sys,[-3 4],"Name","Error");
% Define MFs for input E
sys = addmf(sys,"Error","trimf",[-3 -3 -1], 'Name',"NB");
sys = addmf(sys,"Error","trimf",[-3 -3 -0.25],'Name',"NM");
sys = addmf(sys,"Error","trimf",[-0.5 -0.5 0],'Name',"NS");
sys = addmf(sys,"Error","trimf",[-0.25 0 0.25],'Name',"ZE");
sys = addmf(sys,"Error","trimf",[0 0.5 0.5],'Name',"PS");
sys = addmf(sys,"Error","trimf",[0.25 2 2],'Name',"PM");
sys = addmf(sys,"Error","trimf",[1 4 4],'Name',"PB");

% Define input variable input change-in-error "CE" to fis 
%with range [-1 1]
sys = addInput(sys,[-0.03 0.04],"Name","Change_in_error");
% Define MFs for input CE
sys = addmf(sys,"Change_in_error","trimf",[-0.03 -0.03 -0.02],'Name',"NB");
sys = addmf(sys,"Change_in_error","trimf",[-0.03 -0.02 -0.01],'Name',"NM");
sys = addmf(sys,"Change_in_error","trimf",[-0.02 -0.01 0],'Name',"NS");
sys = addmf(sys,"Change_in_error","trimf",[-0.01 0 0.01],'Name',"ZE");
sys = addmf(sys,"Change_in_error","trimf",[0 0.01 0.02],'Name',"PS");
sys = addmf(sys,"Change_in_error","trimf",[0.01 0.02 0.03],'Name',"PM");
sys = addmf(sys,"Change_in_error","trimf",[0.02 0.04 0.04],'Name',"PB");
%% Outputs
%Define output variable Kpid to fis
% with range [.....]
sys = addOutput(sys,[-24 24],"Name","Kp",'MFType',"constant");

%Define MFs Kp1 -> Kp7 for input Kp
sys = addmf(sys, 'Kp',"constant",-24,'Name',"Ln");
sys = addmf(sys, 'Kp',"constant",-16,'Name',"Lm");
sys = addmf(sys, 'Kp',"constant",-8,'Name',"Ls");
sys = addmf(sys, 'Kp',"constant",0,'Name',"Z");
sys = addmf(sys, 'Kp',"constant",8,'Name',"Ss");
sys = addmf(sys, 'Kp',"constant",16,'Name',"Sm");
sys = addmf(sys, 'Kp',"constant",24,'Name',"Sp");

sys = addOutput(sys,[-3 3],"Name","Ki",'MFType',"constant");
%Define MFs Kd1 -> Kd7 for input Kd
sys = addmf(sys, 'Ki','constant',-3,'Name',"Ln");
sys = addmf(sys, 'Ki','constant',-2,'Name',"Lm");
sys = addmf(sys, 'Ki','constant',-1,'Name',"Ls");
sys = addmf(sys, 'Ki','constant',0,'Name',"Z");
sys = addmf(sys, 'Ki','constant',1,'Name',"Ss");
sys = addmf(sys, 'Ki','constant',2,'Name',"Sm");
sys = addmf(sys, 'Ki','constant',3,'Name',"Sp");


sys = addOutput(sys,[-160 160],"Name","Kd",'MFType',"constant");
%Define MFs Kd1 -> Kd7 for input Kd
sys = addmf(sys, 'Kd','constant',-160,'Name',"Ln");
sys = addmf(sys, 'Kd','constant',-80,'Name',"Lm");
sys = addmf(sys, 'Kd','constant',-40,'Name',"Ls");
sys = addmf(sys, 'Kd','constant',0,'Name',"Z");
sys = addmf(sys, 'Kd','constant',40,'Name',"Ss");
sys = addmf(sys, 'Kd','constant',80,'Name',"Sm");
sys = addmf(sys, 'Kd','constant',160,'Name',"Sp");
%% Add Rule-Base                 
%So the rule-base is interpreted as
% rulelist = [E CE Kp Ki Kd Weight &&=1;||=2]
rulelist = [1 1 2 1 7 1 1;1 2 1 1 7 1 1;1 3 3 1 7 1 1;1 4 4 2 7 1 1;
            1 5 3 1 7 1 1;1 7 2 1 7 1 1;1 7 1 1 7 1 1;
            
            2 1 1 1 7 1 1;2 2 1 1 7 1 1;2 3 1 1 7 1 1;2 4 1 4 7 1 1;
            2 5 2 1 7 1 1;2 6 2 1 6 1 1;2 7 2 2 6 1 1;
            
            3 1 2 2 6 1 1;3 2 2 2 6 1 1;3 3 2 2 6 1 1;3 4 3 2 6 1 1;
            3 5 2 2 4 1 1;3 6 2 2 4 1 1;3 7 3 2 4 1 1;
            
            4 1 7 6 1 1 1;4 2 7 6 2 1 1;4 3 7 7 1 1 1;4 4 7 6 1 1 1;
            4 5 7 6 1 1 1;4 6 7 5 1 1 1;4 7 7 4 1 1 1;
            
            5 1 5 2 3 1 1;5 2 5 2 5 1 1;5 3 4 2 5 1 1;5 4 7 3 1 1 1;
            5 5 7 2 1 1 1;5 6 6 5 2 1 1;5 7 6 2 2 1 1;
            
            6 1 2 4 1 1 1;6 2 2 4 2 1 1;6 3 5 2 2 1 1;6 4 4 5 2 1 1;
            6 5 5 4 1 1 1;6 6 4 3 1 1 1;6 7 4 4 2 1 1;
            
            7 1 7 7 1 1 1;7 2 7 7 1 1 1;7 3 7 6 1 1 1;7 4 7 6 2 1 1;
            7 5 7 7 1 1 1;7 6 7 7 4 1 1;7 7 7 7 2 1 1];
sys = addrule(sys,rulelist);
%%
global t m c k F 
m = 1; k = 5; c = 9;
% Initial states
x = 2; dx = -2;
x0 = [x dx];
X = x0;
%Thoi gian 
tsamp = 0.01;
t = 20;
T = 0;
Runing_time = t/tsamp;
tspan = [0 tsamp];
ts = 0.01;

ref = 0.5; % gia tri mong muon
% he so PID
kp = 12.5; ki = 1.3; kd = 900;

for i = 1:Runing_time  
    % Fix value PID Controller
    time(i) = i*ts;
    tp = i*ts;
    
    if (i==1)
        e(1) = ref;
        erot(1) = 0;
        u(1) = 0;
        upid(1) = 0;
    else
        e(i) = ref - position(i-1);
        erot(i) = e(i) - e(i-1);
        
        u = evalfis([e(i) erot(i)],sys);
        dkp(i) = u(1);
        dki(i) = u(2);
        dkd(i) = u(3);
        upid(i) = (kp+dkp(i))*e(i) + (ki+dki(i))*(e(i) + e(i-1)) + (kd+dkd(i))*(e(i) - e(i-1));
    end
    
    FF(i) = upid(i);
    F = FF(i);
    if (i==1)
        y(1) = 0;
    else
        [t,y] = ode45(@msd,tspan,x0);
        x0 = y(length(y),:);
    end
    T = [T;i*tsamp];
    X = [X;x0];
    position(i) = X(i,1);  
end
%%
for i = 1:Runing_time
    reference(i) = ref;
end
plot(time,reference,'-- b','LineWidth',2);
hold on;
plot(T,X(:,1),'r','LineWidth',2);
xlabel('Time(s)');
ylabel('System output');
ylim([0 5]);
title('Response of the mass position using Fuzzy-PID controller');
%legend('Reference output','FLC');
grid on;
%% Mathematical model function for ODE45's solution
function dxdt = msd(t,z)
    global F m c k
    dxdt_1 = z(2);
    dxdt_2 = (F + 5 - c*z(2) - k*z(1))/m;
    dxdt = [dxdt_1; dxdt_2];
end