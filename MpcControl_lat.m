classdef MpcControl_lat < MpcControlBase
    
    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   x0           - initial state (estimate)
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   u0           - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_segs = ceil(mpc.H/mpc.Ts); % Horizon steps
            N = N_segs + 1;              % Last index in 1-based Matlab indexing

            [nx, nu] = size(mpc.B);
            
            % Targets
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);

            % Initial states
            x0 = sdpvar(nx, 1);
            x0other = sdpvar(nx, 1); % (Ignore this, not used)

            % Input to apply to the system
            u0 = sdpvar(nu, 1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D
            %       are the DISCRETE-TIME MODEL of your system
            %       You can find the linearization steady-state
            %       in mpc.xs and mpc.us.

            A = mpc.A;
            B = mpc.B;
            Q = 1*eye(2);% Adjust
            R = 1 ; % Adjust
            M = [1;-1]; m = [0.5236; 0.5236];
            F = [1 0;-1 0;0 1;0 -1]; f = [3.5; 0.5; 0.0873; 0.0873];
            [K,Qf,~] = dlqr(A,B,Q,R);
            K = -K;

            % Compute maximal invariant set
            Xf = polytope([F;M*K],[f;m]);
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
            figure
            hold on; grid on;
            plot(polytope([F;M*K],[f;m]),'g'); 
            plot(Xf,'r');
            legend('State & Input constraints set', 'Maximal invariant set');
            xlabel('y position [m]'); ylabel('Heading angle [rad]');

            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            x = sdpvar(nx,N,'full');
            u = sdpvar(nu,N,'full');
            obj = 0;
            con = [];
            for i = 1:N-1               
                con = [con,(x(:,i+1) == mpc.f_xs_us + mpc.A*(x(:,i) - mpc.xs) + B*(u(:,i) - mpc.us))]; % System dynamics
                con = [con , (M*(u(:,i)) <= m)]; % Input constraints
                con = [con , (F*(x(:,i)) <= f)]; % Input constraints
                obj = obj + (x(:,i) - x_ref)'*Q*(x(:,i) - x_ref) + (u(:,i) - u_ref)'*R*(u(:,i) - u_ref); % Cost function
            end
            con = [con , (Ff*(x(:,N) - x_ref) <= ff)];
            obj = obj + (x(:,N)-x_ref)'*Qf*(x(:,N) - x_ref); % Terminal weight

            % Replace this line and set u0 to be the input that you
            % want applied to the system. Note that u0 is applied directly
            % to the nonlinear system. You need to take care of any 
            % offsets resulting from the linearization.
            % If you want to use the delta formulation make sure to
            % substract mpc.xs/mpc.us accordingly.
            con = [con , ( u0 == u(:,1))];
            con = [con , ( x0 == x(:,1))];

            % Pass here YALMIP sdpvars which you want to debug. You can
            % then access them when calling your mpc controller like
            % [u, X, U] = mpc_lat.get_u(x0, ref);
            % with debugVars = {X_var, U_var};
            debugVars = {x,u};

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {x0, x_ref, u_ref, x0other}, {u0, debugVars{:}});
        end
        
        % Computes the steady state target which is passed to the
        % controller
        function [xs_ref, us_ref] = compute_steady_state_target(mpc, ref)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            % OUTPUTS
            %   xs_ref, us_ref - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Steady-state system
            A = mpc.A;
            B = mpc.B;

            % Linearization steady-state
            xs = mpc.xs;
            us = mpc.us;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            Xss = sdpvar(2,1);
            Uss = sdpvar(1,1);
            y_ref = sdpvar(1,1);
            M = [1;-1]; m = [0.5236; 0.5236];
            F = [1 0;-1 0;0 1;0 -1]; f = [3.5; 0.5; 0.0873; 0.0873];
            % Constraints
            con = [];
            con = [con, xs + A*(Xss-xs) + B*(Uss-us) == Xss]; % Stationnary state
            con = [con, F*Xss <= f];  % State constraints
            con = [con, M*Uss <= m];  % Input constraints
            con = [con, mpc.C(2,2)*Xss(1) == y_ref];
            obj = Uss^2;
            opt = optimizer(con, obj, sdpsettings('solver', 'gurobi'), y_ref, [Xss; Uss]);
            res = opt(ref);
            xs_ref = double(res(1:2));
            us_ref = double(res(3));
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
