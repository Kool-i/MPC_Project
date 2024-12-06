car = Car(Ts);
Vs = 120/3.6; % 120 km/h
[xs, us] = car.steady_state(Vs); % Compute steady−state for which f s(xs,us) = 0
sys = car.linearize(xs, us); % Linearize the nonlinear model around xs, us
rank(ctrb(sys.A,sys.B))
rank(obsv(sys.A,sys.C))
[sys_lon, sys_lat] = car.decompose(sys);
Ts = 1/10;
[fd_xs_us, Ad, Bd, Cd, Dd] = Car.c2d_with_offset(sys, Ts);