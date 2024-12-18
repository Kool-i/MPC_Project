classdef v10dec_MpcControl_lon < MpcControlBase
    
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
            A = mpc.A;
            B = mpc.B;
            Q = 10 * eye(2); % Adjust
            R = 1 ; % Adjust

            % Constraints
            % u in U = { u | Mu <= m }
            M = [1;-1]; m = [1; 1];
            
            [K,Qf,~] = dlqr(A,B,Q,R);
            % MATLAB defines K as -K, so invert its signal
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
            size_Xf = size(Xf)
            
            % Visualizing the sets
            figure
            hold on; grid on;
            plot(polytope([],[]),'g'); plot(Xf,'r'); % A ajuster
            xlabel('position'); ylabel('velocity'); % A ajuster

            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D
            %       are the DISCRETE-TIME MODEL of your system.
            %       You can find the linearization steady-state
            %       in mpc.xs and mpc.us.
            
            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            x = sdpvar(nx,N,'full');
            u = sdpvar(nu,N,'full');
            obj = 0;
            con = [];

            % Replace this line and set u0 to be the input that you
            % want applied to the system. Note that u0 is applied directly
            % to the nonlinear system. You need to take care of any 
            % offsets resulting from the linearization.
            % If you want to use the delta formulation make sure to
            % substract mpc.xs/mpc.us accordingly.
            con = con + ( u0 == u(:,1));
            con = con + ( x0 == x(:,1));

            % Pass here YALMIP sdpvars which you want to debug. You can
            % then access them when calling your mpc controller like
            % [u, X, U] = mpc_lon.get_u(x0, ref);
            % with debugVars = {X_var, U_var};
            debugVars = {x, u};
            sz = size(x)
            sz = size(u)
            sz = size(mpc.us)
            sz = size(A)
            sz = size(B)
            a = A
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1:N-1
                % con = con + (x(:,i+1) == mpc.f_xs_us(2) + A*(x(:,i) - mpc.xs(2)) + B*(u(:,i) - mpc.us)); % System dynamics
                con = con + (x(:,i+1) == mpc.xs + A*(x(:,i) - mpc.xs(2)) + B*(u(:,i) - mpc.us)); % System dynamics
                con = con + (M*u(:,i) <= m); % Input constraints
                obj = obj + (x(2,i) - V_ref)'*Q(2,2)*(x(2,i)-V_ref) + (u(:,i) - u_ref)'*R*(u(:,i) - u_ref); % Cost function
            end
            con = con + (Ff*(x(:,N)) <= ff);
            obj = obj + (x(2,N)-V_ref)'*Qf(2,2)*(x(2,N) - V_ref); % Terminal weight
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
            Vs_ref = 0;
            us_ref = 0;
            Xs = sdpvar(1,1);
            Us = sdpvar(1,1);
            % u in U = { u | Mu <= m }
            M = [1;-1]; m = [1; 1];
            % Constraints
            con = [A*(Xs-xs) + B*(Us-us) == Xs-xs]; % Stationnary state
            con = [con, M*(Us) <= m];  % Input constraints
            con = [con, mpc.C*(Xs) == ref];
            obj = (Us-us)'*(Us-us);

            % %v10dec
            % obj = 0;
            % con = [xs == 0, us == 0];
        
            % Compute the steady-state target
            % optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {Xs, Us});
            %optimize(con, obj, sdpsettings('solver', 'gurobi'));
            % Sorties
            %Vs_ref = double(Xs);
            %us_ref = double(Us);

            A = mpc.A(2,2);
            B = mpc.B(2,1);
            Vs_ref = ref;
            us_ref=A*(1-Vs_ref)/B;
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
