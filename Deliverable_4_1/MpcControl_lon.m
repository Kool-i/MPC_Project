classdef MpcControl_lon < MpcControlBase
    
    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   x0           - initial state (estimate)
            %   V_ref, u_ref - reference state/input
            %   d_est        - disturbance estimate
            %   x0other      - initial state of other car
            % OUTPUTS
            %   u0           - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_segs = ceil(mpc.H/mpc.Ts); % Horizon steps
            N = N_segs + 1;              % Last index in 1-based Matlab indexing

            [nx, nu] = size(mpc.B);
            
            % Targets
            V_ref = sdpvar(1);
            u_ref = sdpvar(1);

            % Disturbance estimate (Ignore this before Todo 4.1)
            d_est = sdpvar(1);

            % Initial states
            x0 = sdpvar(nx, 1);
            x0other = sdpvar(nx, 1); % (Ignore this before Todo 5.1)

            % Input to apply to the system
            u0 = sdpvar(nu, 1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D
            %       are the DISCRETE-TIME MODEL of your system.
            %       You can find the linearization steady-state
            %       in mpc.xs and mpc.us.

            A = mpc.A(2,2);
            B = mpc.B(2,1);
            Q = 3;
            R = 1;
            M = [1;-1]; m = [1; 1];
            [K,Qf,~] = dlqr(A,B,Q,R);
            K = -K;

            % Compute maximal invariant set
            Xf = polytope([M*K],[m]);
            Acl = [A+B*K];
            while 1
                prevXf = Xf;
                [T,t] = double(Xf);
                preXf = polytope(T*Acl,t);
                Xf = intersect(Xf, preXf);
                if isequal(prevXf, Xf)
                    break
                end
            end
            [Ff,ff] = double(Xf);

            % Visualizing the sets
            % figure
            % hold on; grid on;
            % plot(Xf,'r');
            % legend('Maximal invariant set');
            % xlabel('x position [m]'); ylabel('Heading angle [rad]');
            
            x = sdpvar(1,N,'full');
            u = sdpvar(nu,N,'full');
            obj = 0;
            con = [];
            for i = 1:N-1
                con = [con , (x(1,i+1) == A*(x(1,i)) + B*u(:,i))]; % System dynamics
                con = [con , (M*(u(:,i) + u_ref) <= m)]; % Input constraints
                obj = obj + (x(1,i))'*Q*(x(1,i)) + (u(:,i))'*R*(u(:,i)); % Cost function
            end
            con = [con , (Ff*(x(1,N)) <= ff)];
            obj = obj + (x(1,N))'*Qf*(x(1,N)); % Terminal weight

            % Replace this line and set u0 to be the input that you
            % want applied to the system. Note that u0 is applied directly
            % to the nonlinear system. You need to take care of any 
            % offsets resulting from the linearization.
            % If you want to use the delta formulation make sure to
            % substract mpc.xs/mpc.us accordingly.
            con = [con , ( u(:,1) == u0 - u_ref)];
            con = [con , ( x(1,1) == x0(2) - V_ref)];

            % Pass here YALMIP sdpvars which you want to debug. You can
            % then access them when calling your mpc controller like
            % [u, X, U] = mpc_lon.get_u(x0, ref);
            % with debugVars = {X_var, U_var};
            debugVars = {x,u};

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {x0, V_ref, u_ref, d_est, x0other}, {u0, debugVars{:}});
        end
        
        % Computes the steady state target which is passed to the
        % controller
        function [Vs_ref, us_ref] = compute_steady_state_target(mpc, ref, d_est)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            %   d_est  - disturbance estimate (Ignore before Todo 4.1)
            % OUTPUTS
            %   Vs_ref, us_ref - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Steady-state subsystem
            A = mpc.A(2, 2);
            B = mpc.B(2, 1);

            % Subsystem linearization steady-state
            xs = mpc.xs(2);
            us = mpc.us;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

            r = [-A*xs - B*(us - d_est) + xs; ref];
            M = [1 - A, -B; 1, 0];
            sol = (M^-1)*r;
            Vs_ref = sol(1);
            us_ref = sol(2);

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
