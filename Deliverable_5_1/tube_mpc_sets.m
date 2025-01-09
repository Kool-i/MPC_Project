addpath(genpath("Utils"))
addpath(genpath("Part_5"))

clear
close all

Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us_all] = car.steady_state(120 / 3.6);
sys = car.linearize(xs, us_all);
[sys_lon, ~] = car.decompose(sys);

[Ad, Bd, ~, ~] = ssdata(c2d(sys_lon, Ts));
Bd = -Bd;
us = sys_lon.UserData.us;

% PARAMETERS %
xsafe = 7;
Q = diag([3,1]);
R = 1.0;

[K,Qf,~] = dlqr(Ad,Bd,Q,R);
K = -K;
Acl = [Ad+Bd*K];

WM = [1 0;0 1;-1 0;0 -1]; wm = -[Bd(1)*(us+0.5);Bd(2)*(us+0.5);-Bd(1)*(us-0.5);-Bd(2)*(us-0.5)];
Wf = Polyhedron('A',WM,'b',wm);

% Compute minimum robust invariant set
Epsilon = Polyhedron('A',[1, 0; -1, 0; 0, 1; 0, -1],'b',zeros(4,1));
figure(1)
hold on; grid on;

i=0;
while 1
    Epsilon.plot('alpha', 0, 'linewidth', 0.5);

    Wft = Wf.affineMap(Acl^i);
    Epsilonnew = plus(Epsilon, Wft);
    Epsilonnew = Epsilonnew.minHRep;

    disp(norm(Acl^i));
    i=i+1;
    if norm(Acl^i) < 0.01 || i == 100
        break
    end
    Epsilon = Epsilonnew;
end
Epsilon.plot('alpha', 0, 'edgecolor', 'r', 'linewidth', 2);
legend('Epsilon');
xlabel('delta-x position [m]'); ylabel('delta-V speed [m/s]');

F = [-1 0]; f = [-6+xsafe];
M = [1;-1]; m = [1; 1];

X = Polyhedron('A',F,'b',f);
U = Polyhedron('A',M,'b',m);
Xtilde = minus(X, Epsilon);
KEpsilon = Epsilon.affineMap(K);
Utilde = minus(U, KEpsilon);

%Visualizing the sets
figure(2)
hold on; grid on;
Utilde.plot('alpha', 0, 'edgecolor', 'b', 'linewidth', 2);
legend('Utilde');
xlabel('uT');

figure(3)
hold on; grid on;
Xtilde.plot('alpha', 0, 'edgecolor', 'b', 'linewidth', 2);
X.plot('alpha', 0, 'edgecolor', 'r', 'linewidth', 2);
legend('Xtilde','X', location="southeast");
xlabel('delta-x position [m]'); ylabel('delta-V speed [m/s]');


save("tube_mpc_data.mat", 'Epsilon', 'Xtilde', 'Utilde', 'K', 'Qf', 'Q', 'R', 'xsafe');
