classdef NmpcControl_overtake < handle

    properties
        % The NMPC problem
        opti

        % Problem parameters
        x0, ref, x0other

        % Most recent problem solution
        sol

        % The input that you want to apply to the system
        u0

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add any variables you would like to read to debug here
        % and then store them in the NmpcControl function below.
        % e.g., you could place X here and then add obj.X = X
        % in the NmpcControl function below.
        % 
        % After solving the problem, you can then read these variables 
        % to debug via
        %   nmpc.sol.value(nmpc.X)
        % 
        X, U, Xother
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end

    methods
        function obj = NmpcControl_overtake(car, H)

            import casadi.*

            N_segs = ceil(H/car.Ts); % Horizon steps
            N = N_segs + 1;          % Last index in 1-based Matlab indexing

            nx = 4;
            nu = 2;

            % Define the NMPC optimization problem
            opti = casadi.Opti();
            
            % Parameters (symbolic)
            obj.x0 = opti.parameter(nx, 1);       % initial state
            obj.ref = opti.parameter(2, 1);       % target y, velocity
            obj.x0other = opti.parameter(nx, 1);  % initial state of other car

            % SET THIS VALUE TO BE YOUR CONTROL INPUT
            obj.u0 = opti.variable(nu, 1);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

            % Define your problem using the opti object created above


            obj.X = opti.variable(nx, N); 
            obj.Xother = opti.variable(nx, N);
            % Xother(:,1) = obj.x0other;
            obj.U = opti.variable(nu, N-1);
            Matrix_H = zeros(2,2);
            Matrix_H(1,1) = 1/(9*9);
            Matrix_H(2,2) = 1/(3*3);

            % change this line accordingly
            opti.subject_to(obj.U(:, 1) == obj.u0);
            opti.subject_to(obj.X(:,1) == obj.x0)
            opti.subject_to(obj.Xother(:,1) == obj.x0other)

            % Cost function
            Q = diag([0, 0.001, 2, 150]);    % State tracking cost [x, y, θ, V]
            R = diag([2, 1]);           % Control effort cost [δ, uT]
            cost = 0;
            
            for k = 1:N-1
                % State cost (tracking reference)
                % x_err = obj.X(:, k) - [0; obj.ref(1); 0; obj.ref(2)];
                x_err = obj.X(:, k) - [0; 0; 0; obj.ref(2)];
                cost = cost + x_err' * Q * x_err + obj.U(:, k)' * R * obj.U(:, k);
                opti.subject_to(obj.Xother(:, k+1) == [obj.Xother(1, k) + car.Ts * obj.Xother(4, k);
                                       obj.Xother(2, k);
                                       obj.Xother(3, k);
                                       obj.Xother(4, k)]);

            end
            opti.minimize(cost);

            % Constraints
            % Initial condition
            opti.subject_to(obj.X(:, 1) == obj.x0);

            % System dynamics constraints using RK4 integration
            for k = 1:N-1
                x_next = obj.rk4_integrate(car, obj.X(:, k), obj.U(:, k), car.Ts);
                opti.subject_to(obj.X(:, k+1) == x_next);
                opti.subject_to((obj.X(1:2, k+1)-obj.Xother(1:2, k+1))'*Matrix_H*(obj.X(1:2, k+1)-obj.Xother(1:2, k+1)) >= 1); % Heading constraint
            end

            % State constraints
            opti.subject_to(-0.5 <= obj.X(2, :) <= 3.5); % Lane constraints
            opti.subject_to(-0.0873 <= obj.X(3, :) <= 0.0873); % Heading constraint
            

            % Input constraints
            opti.subject_to(-0.5236 <= obj.U(1, :) <= 0.5236); % Steering limit
            opti.subject_to(-1.0 <= obj.U(2, :) <= 1.0); % Throttle limit 


            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Store the defined problem to solve in get_u
            obj.opti = opti;

            % Setup solver
            options = struct;
            options.ipopt.print_level = 0;
            options.print_time = 0;
            options.expand = true;
            obj.opti.solver('ipopt', options);
        end

        function x_next = rk4_integrate(obj, car, x, u, dt)
            % Runge-Kutta 4th order integration for system dynamics
            k1 = car.f(x, u);
            k2 = car.f(x + dt/2 * k1, u);
            k3 = car.f(x + dt/2 * k2, u);
            k4 = car.f(x + dt * k3, u);

            % Compute next state
            x_next = x + dt / 6 * (k1 + 2*k2 + 2*k3 + k4);
        end

        function u = get_u(obj, x0, ref, x0other)

            if nargin < 4
                x0other = zeros(4, 1);
            end

            % Compute solution from x0
            obj.solve(x0(1:4), ref, x0other(1:4));
            u = obj.sol.value(obj.u0);   
        end

        function solve(obj, x0, ref, x0other)

            % Pass parameter values
            obj.opti.set_value(obj.x0, x0);
            obj.opti.set_value(obj.ref, ref);
            obj.opti.set_value(obj.x0other, x0other);

            obj.sol = obj.opti.solve();   % actual solve
            
            % Set warm start for next solve
            obj.opti.set_initial(obj.sol.value_variables());
            obj.opti.set_initial(obj.opti.lam_g, obj.sol.value(obj.opti.lam_g));
        end
    end
end
