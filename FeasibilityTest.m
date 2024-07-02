close all; 
clear all; 
clc; 
 
% Linearized CTS System: 
% Assuming a Rocket length of 1 meter, so L = 0.5 meters, Iz = 54 Kg*m^2   
% Parameters at point to linearize about
Thrust = 1000; %N 
V_e = 10000; %m/s
Iz = 54;  %Kg*m^2
g = 9.8; %m/s^2
mass = 100; %kg
L = 0.5; %m
Ac = [1000/(mass*V_e), 0, -1000/mass, 0, 0;  
    0, 1000/(mass*V_e), 0, 0, -1000/(mass^2);  
    0, 0, 0, 1, 0;  
    0, 0, 0, 0, 0;  
    0, 0, 0, 0, 0]; 
Bc = [0, 1000/mass;  
    1/mass, 0;  
    0, 0;  
    0, (L*1000)/Iz;  
    -1/V_e, 0]; 
% Euler's Method
T = 0.1;
A = eye(5) + T .* Ac;
B = T .* Bc;
% Q is chosen as cost of x state matrix, R is chosen as cost of u input matrix 
Q = eye(5); 
R = eye(2); 
%P_inf is the solution to the discrete algebraic Riccati equation 
[P_inf, K, L] = idare(A, B, Q, R, [], []); 
P = P_inf; 
% D is our affine term, system in the form A*x + B*u + D 
Dc = [0; 1000/mass - g; 0; 0; -1000/V_e]; 
D = T.*Dc;  
% X state polyhedron, lower/upper bounds for X states, Vx, Vy, Theta, There_dot, m 
X = Polyhedron('lb', [-100; -100; -pi/6; -pi/3; 10], 'ub', [100; 100; pi/6; pi/3; 400]); 
% U state polyhedron, lower/upper bounds for U states, Thrust T and Delta 
U = Polyhedron('lb', [-1300; -pi/3], 'ub', [1300; pi/3]); 
%Finding our CI set and labeling it C for our LTI system so Xf = C 
system = LTISystem('A', A, 'B', B, 'f', D); 
system.x.min = [-10; -10; -pi/6; -pi/3; 0]; 
system.x.max = [10; 10; pi/6; pi/3; 200]; 
system.u.min = [-1300; -pi/3]; 
system.u.max = [1300; pi/3]; 
InvSet = system.invariantSet(); 
C = InvSet; 
% Initial condition goes here
x0_list = [0; 0; 0; 0; 400]; 
% T timestep for plotting
T = 200; 
% Arrays to plot 
x_traj = zeros(5, T+1); 
u_traj = zeros(2, T); 
%For each horizon 
for N = [2, 7] 
        %trajectory of x starts with x0
        x_traj(:,1) = x0_list(:, 1); 
        for t = 1:T 
            %batch method to get Abatch, Bbatch, Qbatch, Rbatch 
            [Abatch, Bbatch, Qbatch, Rbatch] = BatchMatrices(A, B, Q, R, P, N); 
            %we want to construct Dbatch to be similar to Bbatch except with 
            %I instead of B, so it looks like 
            % [0 0      ...0] * [D] 
            % [I 0      0..0]   [.] 
            % [A I      0..0]   [.]   
            % [A^(N-1) ... I]   [D] 
            [~, Dbatch, ~, ~] = BatchMatrices(A, eye(5), Q, R, P, N); 
            Dbatch = Dbatch*kron(ones(N,1), D); 
            %H from batch method definition 
            H = 2 * (Bbatch' * Qbatch * Bbatch + Rbatch); 
            %f from batch method definition, with affine term it will be: 
            f = 2 * (Bbatch' * Qbatch' * (Abatch * x_traj(:, t) + Dbatch)); 
            %HBaru is Kron I NxN Hu 
            Hu = kron(eye(N), U.A); 
            %Kbaru is Kron I Nx1 Ku 
            Ku = kron(ones(N,1), U.b); 
            %HBarx is Kron I NxN Hx and then Hf at the end 
            Hx = blkdiag(kron(eye(N), X.A), C.A); 
            %Kbarx is kron I Nx1 kx and then kf at the end 
            Kx = [kron(ones(N,1), X.b); C.b]; 
            % quad prog to get optimal u value 
            [u_opt, ~] = quadprog(H, f, [Hx * Bbatch; Hu], [Kx-Hx * Abatch * x_traj(:, t)-Hx * Dbatch; Ku], [], []); 
            %if outside feesible region 
            if isempty(u_opt) 
                disp('State or input went out of bounds.'); 
            else 
                %update u from the quadprog result (optimal solution) 
                u_traj(:,t) = [u_opt(1); u_opt(2)]; 
                %update X using Ax + Bu + D 
                x_next = A * x_traj(:, t) + B * u_traj(:,t) + D; 
                %update x trajectory to the next x 
                x_traj(:, t+1) = x_next;   
            end 
        end 
        %Plot 
        figure; 
        subplot(7, 1, 1); 
        plot(0:T, x_traj(1,:), 'b', 'LineWidth', 2); 
        xlabel('Time'); 
        ylabel('X1 (Vx)');
        hold on; 
        title(['N: [',num2str(N), ']']); 
        subplot(7, 1, 2); 
        plot(0:T, x_traj(2,:), 'b', 'LineWidth', 2); 
        xlabel('Time'); 
        ylabel('X2 (Vy)'); 
        subplot(7, 1, 3); 
        plot(0:T, x_traj(3,:), 'b', 'LineWidth', 2); 
        xlabel('Time');
        ylabel('X3 (Theta)'); 
        subplot(7, 1, 4); 
        plot(0:T, x_traj(4,:), 'b', 'LineWidth', 2); 
        xlabel('Time'); 
        ylabel('X4 (ThetaDot)'); 
        subplot(7, 1, 5); 
        plot(0:T, x_traj(5,:), 'b', 'LineWidth', 2); 
        xlabel('Time'); 
        ylabel('X5 (mass)'); 
        subplot(7, 1, 6); 
        plot(0:T-1, u_traj(1,:), 'r', 'LineWidth', 2); 
        xlabel('Time'); 
        ylabel('U1 (Thrust)'); 
        subplot(7, 1, 7); 
        plot(0:T-1, u_traj(2,:), 'r', 'LineWidth', 2); 
        xlabel('Time'); 
        ylabel('U2 (Delta)'); 
        grid on; 
    end   
 




% Define number of random initial conditions
num_random_conditions = 1500;

% Define state variable bounds
x_min = [-10; -10; -pi/6; -pi/3; 10]; % Lower bounds for states
x_max = [10; 10; pi/6; pi/3; 200]; % Upper bounds for states

% Ensure 'T' is an integer
T = round(15); % Ensure integer value

% Store results for tracking if they are within feasible regions
out_of_bounds_count = 0;

% Loop to generate and simulate for each random initial condition
for i = 1:num_random_conditions
    % Generate random initial condition
    x0 = x_min + (x_max - x_min) .* rand(5, 1); % Randomized initial condition
    
    % Reset trajectory arrays
    x_traj = zeros(5, T + 1); % Ensure integer array size
    u_traj = zeros(2, T); % Ensure integer array size
    x_traj(:, 1) = x0; 
    
    % Loop through time steps
    for t = 1:T
        % Calculate batch matrices and optimal control inputs
        [Abatch, Bbatch, Qbatch, Rbatch] = BatchMatrices(A, B, Q, R, P, N); 
        [~, Dbatch, ~, ~] = BatchMatrices(A, eye(5), Q, R, P, N);
        Dbatch = Dbatch * kron(ones(N,1), D);
        
        % Calculate H and f for quadratic programming
        H = 2 * (Bbatch' * Qbatch * Bbatch + Rbatch); 
        f = 2 * (Bbatch' * Qbatch' * (Abatch * x_traj(:, t) + Dbatch));
        
        % Define constraints
        Hu = kron(eye(N), U.A); 
        Ku = kron(ones(N,1), U.b);
        Hx = blkdiag(kron(eye(N), X.A), C.A);
        Kx = [kron(ones(N,1), X.b); C.b];
        
        % Solve quadratic programming
        [u_opt, ~] = quadprog(H, f, [Hx * Bbatch; Hu], [Kx-Hx * Abatch * x_traj(:, t)-Hx * Dbatch; Ku], [], []);
        
        % Check if control inputs are feasible
        if isempty(u_opt)
            out_of_bounds_count = out_of_bounds_count + 1; % Out-of-bounds trajectory
            break;
        else
            % Update control and state trajectories
            u_traj(:, t) = [u_opt(1); u_opt(2)];
            x_next = A * x_traj(:, t) + B * u_traj(:, t) + D;
            x_traj(:, t + 1) = x_next; % Ensure valid indexing
        end
    end
    
    % Plot control and state trajectories if within feasible bounds
    % if t == T % Only plot if all steps were within feasible region
    %     figure; 
    %     subplot(7, 1, 1); 
    %     plot(0:T, x_traj(1,:), 'b', 'LineWidth', 2); 
    %     xlabel('Time'); 
    %     ylabel('X1 (Vx)'); 
    %     subplot(7, 1, 2); 
    %     plot(0:T, x_traj(2,:), 'b', 'LineWidth', 2); 
    %     xlabel('Time'); 
    %     ylabel('X2 (Vy)'); 
    %     subplot(7, 1, 3); 
    %     plot(0:T, x_traj(3,:), 'b', 'LineWidth', 2); 
    %     xlabel('Time'); 
    %     ylabel('X3 (Theta)'); 
    %     subplot(7, 1, 4); 
    %     plot(0:T, x_traj(4,:), 'b', 'LineWidth', 2); 
    %     xlabel('Time'); 
    %     ylabel('X4 (ThetaDot)'); 
    %     subplot(7, 1, 5); 
    %     plot(0:T, x_traj(5,:), 'b', 'LineWidth', 2); 
    %     xlabel('Time'); 
    %     ylabel('X5 (mass)'); 
    %     subplot(7, 1, 6); 
    %     plot(0:T - 1, u_traj(1,:), 'r', 'LineWidth', 2); 
    %     xlabel('Time'); 
    %     ylabel('U1 (Thrust)'); 
    %     subplot(7, 1, 7); 
    %     plot(0:T - 1, u_traj(2,:), 'r', 'LineWidth', 2); 
    %     xlabel('Time'); 
    %     ylabel('U2 (Delta)'); 
    %     grid on;
    % end
end

fprintf('Number of trajectories out of bounds: %d in %d\n', out_of_bounds_count, num_random_conditions);
