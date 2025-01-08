addpath(genpath("Utils"))
addpath(genpath("Part_5"))

clear
close all

Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(120 / 3.6);
sys = car.linearize(xs, us);
[sys_lon, sys_lat] = car.decompose(sys);

% Design MPC controller
H_lon = 10; % Horizon length in seconds
H_lat = 3; % Horizon length in seconds
mpc_lon = MpcControl_lon(sys_lon, Ts, H_lon);
mpc_lat = MpcControl_lat(sys_lat, Ts, H_lat);
mpc = car.merge_lin_controllers(mpc_lon, mpc_lat);

estimator = LonEstimator(sys_lon, Ts);

params = {};
params.Tf = 25;
params.myCar.model = car;
params.myCar.x0 = [0 0 0 100/3.6]';
params.myCar.u = @mpc.get_u;
params.myCar.ref = [0 120/3.6]';
params.otherCar.model = car;
params.otherCar.x0 = [15 0 0 100/3.6]';
params.otherCar.u = car.u_const(100/3.6);
result = simulate(params);
visualization(car, result);