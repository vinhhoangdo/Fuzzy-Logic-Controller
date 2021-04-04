%% Mass - Spring - Damper using FLC 
clear all, close all, clc
%% New FIS GUI
%Create new fis
sys = newfis("MassSpringDamper");
%% Inputs
% Define input variable error "E" to fis 
% with range [-4 4]
sys = addvar(sys,'input', 'E', [-1.5 0.2]);
% Define MFs for input E
sys = addmf(sys, 'input',1,'NB','trimf',[-1.5 -1.5 0]);
% sys = addmf(sys, 'input',1,'NS','trimf',[-4 -2 0]);
sys = addmf(sys, 'input',1,'ZE','trimf',[-0.004 0 0.002]);
% sys = addmf(sys, 'input',1,'PS','trimf',[0 2 4]);
sys = addmf(sys, 'input',1,'PB','trimf',[0 0.2 0.2]);

% Define input variable input change-in-error "CE" to fis 
% with range [-1 1]
sys = addvar(sys, 'input','CE',[-0.1 0.2]);
% Define MFs for input CE
sys = addmf(sys, 'input',2,'NB','trimf',[-0.1 -0.1 0]);
% sys = addmf(sys, 'input',2,'NS','trimf',[-0.5 -0.25 0]);
sys = addmf(sys, 'input',2,'ZE','trimf',[-0.0001 0 0.001]);
% sys = addmf(sys, 'input',2,'PS','trimf',[0 0.25 0.5]);
sys = addmf(sys, 'input',2,'PB','trimf',[0 0.2 0.2]);
%% Outputs
% Define output variable U to fis
% with range [-70 15]
sys = addvar(sys,'output','U',[0 7]);
% Define MFs for input U
sys = addmf(sys, 'output',1,'NB','trimf',[0 0 3]);
sys = addmf(sys, 'output',1,'NS','trimf',[2.5 3 4.5]);
sys = addmf(sys, 'output',1,'ZE','trimf',[3 4 5]);
sys = addmf(sys, 'output',1,'PS','trimf',[4 6 6]);
sys = addmf(sys, 'output',1,'PB','trimf',[5 7 7]);

%% Add Rule-Base                 
% So the rule-base is interpreted as
% rulelist = [E CE U Weight &&=1;||=2)]
rulelist = [1 1 1 1 1;1 1 2 1 1;1 2 3 1 1;1 2 4 1 1; 1 3 3 1 1; 1 3 4 1 1;
            2 1 1 1 1; 2 1 2 1 1;2 2 2 1 1;2 2 3 1 1;2 2 4 1 1;2 2 5 1 1;2 3 4 1 1; 2 3 5 1 1;
            3 1 1 1 1;3 2 3 1 1; 3 2 4 1 1; 3 3 3 1 1];
        
sys = addrule(sys,rulelist);
%% Initial parameters & Reference
global t m c k F 
m = 1; c = 5; k = 9;
% Initial states
x = 2; dx = -2;
x0 = [x dx];
X = x0;
%Sample time
tsamp = 0.01;
t = 20;
T = 0;
Runing_time = t/tsamp; 
tspan = [0 tsamp];
ts = 0.01;
ref = 0.5; %reference

for i = 1:Runing_time  
    % Fuzzy Controller
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
        u(i) = evalfis([e(i) erot(i)],sys);
        %upid(i) = kp*e(i) + ki*(e(i) + e(i-1)) + kd*(e(i) - e(i-1));
    end
    
    FF(i) = u(i);
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
%% Plot step responce
for i = 1:Runing_time
    reference(i) = ref;
end
plot(time,reference,'-- b','LineWidth',2);
hold on;
plot(T,X(:,1),'r','LineWidth',2);
xlabel('Time(s)');
ylabel('System output');
ylim([0 2]);
title('Response ofposition using Fuzzy controller');
legend('Reference output','FLC');
grid on;
%% View MFs and all of rule
% figure
% plotmf(sys,'input',1)
% figure
% plotmf(sys,'input',2)
% figure
% plotmf(sys,'output',1)
% ruleview(sys)
%% Mathematical model function for ODE45's solution
function dxdt = msd(t,z)
    global F m c k
    dxdt_1 = z(2);
    dxdt_2 = (F - c*z(2) - k*z(1))/m;
    dxdt = [dxdt_1; dxdt_2];
end