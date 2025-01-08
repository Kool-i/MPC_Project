addpath(genpath("Utils"))


%% Part 6.1 open-loop
Ts = 1/10; % Sample time
car = Car(Ts); % Ensure car object is correctly created

H = 7; % Horizon in seconds (must be a scalar)
disp(car.Ts); % Debugging check (should print 0.1)

mpc = NmpcControl(car, H); % Ensure both car and H are passed correctly

x0 = [0 0 0 80/3.6]'; % Initial state
ref = [3 100/3.6]'; % Reference state

u = mpc.get_u(x0, ref); % Check open-loop prediction % check if the open−loop prediction is reasonable

%% Part 6.1 lane change and acceleration 80 km/h to 120 km/h
params = {};
params.Tf = 15;
params.myCar.model = car;
params.myCar.x0 = [0 0 0 80/3.6]';
params.myCar.u = @mpc.get_u;
ref1 = [0 80/3.6]';
ref2 = [3 100/3.6]';
params.myCar.ref = car.ref_step(ref1, ref2, 2);
result = simulate(params);
visualization(car, result)


%% Part 6.2 

H = 10;
mpc = NmpcControl_overtake(car, H);
x0_ego = [0 0 0 80/3.6]';
x0_other = [20 0 0 80/3.6]';
ref1 = [0 80/3.6]';
ref2 = [0 100/3.6]';
params = {};
params.Tf = 15;
params.myCar.model = car;
params.myCar.x0 = x0_ego;
params.myCar.u = @mpc.get_u;
params.myCar.ref = car.ref_step(ref1, ref2, 1);
params.otherCar.model = car;
params.otherCar.x0 = x0_other;
params.otherCar.u = car.u_const(80/3.6);
result = simulate(params);
visualization(car, result);

