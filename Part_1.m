%% Todo 1.1

Ts = 1/10;
car = Car(Ts);
% u_T = 0.1;
% delta = 0;
% x = 0;
% y = 0;
% theta = 0;
% V = 0;
u = [delta, u_T]'; % (Assign appropriately)
x = [x, y, theta, V]'; % (Assign appropriately)
x_dot = car.f(x, u)

%% Todo 1.2 |

car = Car(Ts);
Tf = 2.0; % Simulation end time
x0 = [0, 0, deg2rad(-2), 20/3.6]'; % (x, y, theta, V) Initial state
u = [deg2rad(1), 0.7]'; % (delta, u T) Constant input
x = x0 ; % (Assign appropriately)
x_dot = car.f(x, u)
params = {}; % Setup simulation parameter struct
params.Tf = Tf;
params.myCar.model = car;
params.myCar.x0 = x0;
params.myCar.u = u;
result = simulate(params); % Simulate nonlinear model
visualization(car, result);
result.T % Time at every simulation step
result.myCar.X % State trajectory
result.myCar.U % Input trajectory
