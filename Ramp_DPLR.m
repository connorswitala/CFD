%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Perfectly Working, do not change! 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


run('RampGridGenerator.m');

%% Initialize U

M = 2.5;
P = 10000;
T = 300;
gamma = 1.4;
R = 287;
a = sqrt(gamma*R*T);
rho = P/(R*T);
u = M*a;
v = 0;
CFL = 1;

%% Initialize U

W_inlet = [rho; u; v; P];
U_inlet = primToCons(W_inlet);
U_old = zeros(4,Nx,Ny);

for j = 1:Ny
    for i = 1:Nx
        U_old(:,i,j) = U_inlet;
    end
end


%% Start physics solving

% load('ramp_simulation_data.mat', 'U', 'dU', 't', 'n', 'outer_residual')      
% U_old = U;
% t_old = t(n-1);
% 
dU = zeros(4,Nx,Ny); % Initial guess for dU
outer_residual = 10;
n = 1;
t_old = 0;

U = U_old;
dt = Calculate_dt(U, S, Nx, Ny, CFL);

while outer_residual >= 10e-6
    
    dU_old = dU;
    U_old = U;
    dU = SolveOneTimestep(dU_old, U_old, U_inlet, Nx, Ny, S, V, Face_Normal, dt);
    U = U_old + dU;    
    t(n) = t_old + dt;
    t_old = t(n);
    outer_residual(n) = Calculate_Residual(U,V,S,Face_Normal,Nx,Ny);
   

    if abs(mod(n, 1)) == 0
        run('Postprocessingramp.m');
        % save('ramp_simulation_data.mat', 'U', 'dU', 't', 'n', 'outer_residual')      
        look = sprintf('Outer Residual = %e,  dt = %e, Iteration number = %d.\n', outer_residual(n), dt, n);
        disp(look);
    end
    
    dt = Calculate_dt(U, S, Nx, Ny, CFL);
    n = n+1;

end

surface_x_400x200 = squeeze(xcenter(:,1));
save('ramp_simulation_data.mat', 'U', 'dU', 't', 'n', 'outer_residual', 'surface_x_400x200');

%% Plot final state 
W = zeros(4,Nx,Ny);

for i = 1:Nx
    for j = 1:Ny
        [W(:,i,j)] = consToPrim(U(:,i,j));
    end
end

run('Postprocessingramp.m');

%% Function that goes through each lines and calls the propper DPLR for it. Repeats until convergence.

function [dU] = SolveOneTimestep(dU_old, U, U_inlet, Nx, Ny, S, V, Face_Normal, dt)
    
    global_inner_residual = 1;
    dU_new = zeros(4,Nx,Ny);
  
    while global_inner_residual >= 10e-7

        global_inner_residual = 0;

        for i = 1:Nx
            
            if i == 1

                dU_new(:,1,:) = Left_Boundary_Tridiagonal_Solver(dU_old, U, U_inlet, Ny, S, V, Face_Normal, dt);  

            elseif i > 1 && i < Nx

                dU_new(:,i,:) = Inner_Tridiagonal_Solver(dU_old, U, Ny, S, V, Face_Normal, dt, i);
                global_inner_residual = global_inner_residual + (calculateLocalResidual(dU_new, dU_old, U, V, dt, Ny, S, Face_Normal, i))^2;

            elseif i == Nx

                dU_new(:,Nx,:) = Right_Boundary_Tridiagonal_Solver(dU_old, U, Ny, Nx, S, V, Face_Normal, dt);               

            end

        end

        global_inner_residual = sqrt(global_inner_residual);
        dU_old = dU_new;

    end

    dU = dU_old;

end

%% Function that relaxes the left-most line
function[dU_new] = Left_Boundary_Tridiagonal_Solver(dU_old, U, U_inlet, Ny, S, V, Face_Normal, dt)
    
    g = zeros(4,4,Ny); 
    v = zeros(4,1,Ny);
    dU_new = zeros(4,Ny);
    Ei = zeros(4,4);
    Ugi = U_inlet;

    [Ugtw, Ewt] = Wallghost(U(:,1,Ny), Face_Normal(1,4,1,Ny), Face_Normal(2,4,1,Ny));
    
    % Top

    a =  V(1,Ny)/dt*eye(4) ...
                     + (AP(U(:,1,Ny), Ugi, Face_Normal(1,1,1,Ny), Face_Normal(2,1,1,Ny)) + Ei*AM(U(:,1,Ny), Ugi, Face_Normal(1,1,1,Ny), Face_Normal(2,1,1,Ny)))*S(1,1,Ny) ...                     % Left (Inlet BC)
                     + AP(U(:,1,Ny), U(:,2,Ny), Face_Normal(1,2,1,Ny), Face_Normal(2,2,1,Ny))*S(2,1,Ny) ...                                                 % Right
                     + AP(U(:,1,Ny), U(:,1,Ny-1), Face_Normal(1,3,1,Ny), Face_Normal(2,3,1,Ny))*S(3,1,Ny) ...                                              % Bottom
                     + (AP(U(:,1,Ny), Ugtw, Face_Normal(1,4,1,Ny), Face_Normal(2,4,1,Ny)) + Ewt*AM(U(:,1,Ny), Ugtw, Face_Normal(1,4,1,Ny), Face_Normal(2,4,1,Ny)))*S(4,1,Ny);                       % Top (Wall BC) 
    
    c = AM(U(:,1,Ny), U(:,1,Ny-1), Face_Normal(1,3,1,Ny), Face_Normal(2,3,1,Ny))*S(3,1,Ny);
    
    f =  -( AP(U(:,1,Ny), Ugi, Face_Normal(1,1,1,Ny), Face_Normal(2,1,1,Ny))*S(1,1,Ny)*U(:,1,Ny) + AM(U(:,1,Ny), Ugi, Face_Normal(1,1,1,Ny), Face_Normal(2,1,1,Ny))*S(1,1,Ny)*Ugi...                            % Left (Inlet BC at time n)
            + AP(U(:,1,Ny), U(:,2,Ny), Face_Normal(1,2,1,Ny), Face_Normal(2,2,1,Ny))*S(2,1,Ny)*U(:,1,Ny) + AM(U(:,1,Ny), U(:,2,Ny), Face_Normal(1,2,1,Ny), Face_Normal(2,2,1,Ny))*S(2,1,Ny)*U(:,2,Ny)...          % Right (time n)
            + AP(U(:,1,Ny), U(:,1,Ny-1), Face_Normal(1,3,1,Ny), Face_Normal(2,3,1,Ny))*S(3,1,Ny)*U(:,1,Ny) + AM(U(:,1,Ny), U(:,1,Ny-1), Face_Normal(1,3,1,Ny), Face_Normal(2,3,1,Ny))*S(3,1,Ny)*U(:,1,Ny-1)...  % Bottom (time n)
            + AP(U(:,1,Ny), Ugtw, Face_Normal(1,4,1,Ny), Face_Normal(2,4,1,Ny))*S(4,1,Ny)*U(:,1,Ny) + AM(U(:,1,Ny), Ugtw, Face_Normal(1,4,1,Ny), Face_Normal(2,4,1,Ny))*S(4,1,Ny)*Ugtw )...                       % Top (Wall BC at time n)    
            - AM(U(:,1,Ny), U(:,2,Ny), Face_Normal(1,2,1,Ny), Face_Normal(2,2,1,Ny))*S(2,1,Ny)*dU_old(:,2,Ny);                                                     % Right delta U at time n
    
    alpha = a;
    g(:,:,Ny) = alpha\c;
    v(:,1,Ny) = alpha\f;

    % Middle

    for j = (Ny-1):-1:2

        b = AM(U(:,1,j), U(:,1,j+1), Face_Normal(1,4,1,j), Face_Normal(2,4,1,j))*S(4,1,j);
    
        a =  V(1,j)/dt*eye(4) ...
                     + (AP(U(:,1,j), Ugi, Face_Normal(1,1,1,j), Face_Normal(2,1,1,j)) + Ei*AM(U(:,1,j), Ugi, Face_Normal(1,1,1,j), Face_Normal(2,1,1,j)))*S(1,1,j) ...                     % Left (Inlet BC)
                     + AP(U(:,1,j), U(:,2,j), Face_Normal(1,2,1,j), Face_Normal(2,2,1,j))*S(2,1,j) ...                                                 % Right
                     + AP(U(:,1,j), U(:,1,j-1), Face_Normal(1,3,1,j), Face_Normal(2,3,1,j))*S(3,1,j) ...                                              % Bottom
                     + AP(U(:,1,j), U(:,1,j+1), Face_Normal(1,4,1,j), Face_Normal(2,4,1,j))*S(4,1,j);                                                  % Top 
    
        c = AM(U(:,1,j), U(:,1,j-1), Face_Normal(1,3,1,j), Face_Normal(2,3,1,j))*S(3,1,j); 
    
        f =  -( AP(U(:,1,j), Ugi, Face_Normal(1,1,1,j), Face_Normal(2,1,1,j))*S(1,1,j)*U(:,1,j) + AM(U(:,1,j), Ugi, Face_Normal(1,1,1,j), Face_Normal(2,1,1,j))*S(1,1,j)*Ugi...                           % Left (Inlet BC at time n)
                + AP(U(:,1,j), U(:,2,j), Face_Normal(1,2,1,j), Face_Normal(2,2,1,j))*S(2,1,j)*U(:,1,j) + AM(U(:,1,j), U(:,2,j), Face_Normal(1,2,1,j), Face_Normal(2,2,1,j))*S(2,1,j)*U(:,2,j)...            % Right (time n)
                + AP(U(:,1,j), U(:,1,j-1), Face_Normal(1,3,1,j), Face_Normal(2,3,1,j))*S(3,1,j)*U(:,1,j) + AM(U(:,1,j), U(:,1,j-1), Face_Normal(1,3,1,j), Face_Normal(2,3,1,j))*S(3,1,j)*U(:,1,j-1)...    % Bottom (time n)
                + AP(U(:,1,j), U(:,1,j+1), Face_Normal(1,4,1,j), Face_Normal(2,4,1,j))*S(4,1,j)*U(:,1,j) + AM(U(:,1,j), U(:,1,j+1), Face_Normal(1,4,1,j), Face_Normal(2,4,1,j))*S(4,1,j)*U(:,1,j+1) )...    % Top (time n)               
                - AM(U(:,1,j), U(:,2,j), Face_Normal(1,2,1,j), Face_Normal(2,2,1,j))*S(2,1,j)*dU_old(:,2,j);                                                    % Right delta U at time n-1
        
        alpha = a - b*g(:,:,j+1);
        g(:,:,j) = alpha\c;
        v(:,1,j) = alpha\(f-b*v(:,1,j+1));

    end

    % Bottom
    
    [Ugbw, Ewb] = Wallghost(U(:,1,1), Face_Normal(1,3,1,1), Face_Normal(2,3,1,1));
    
    b = AM(U(:,1,1), U(:,1,2), Face_Normal(1,4,1,1), Face_Normal(2,4,1,1))*S(4,1,Ny); 
    
    a =  V(1,1)/dt*eye(4) ... 
                     + (AP(U(:,1,1), Ugi, Face_Normal(1,1,1,1), Face_Normal(2,1,1,1)) + Ei*AM(U(:,1,1), Ugi, Face_Normal(1,1,1,1), Face_Normal(2,1,1,1)))*S(1,1,1) ...                          % Left (Inlet BC)
                     + AP(U(:,1,1), U(:,2,1), Face_Normal(1,2,1,1), Face_Normal(2,2,1,1))*S(2,1,1) ...                                                      % Right
                     + (AP(U(:,1,1), Ugbw, Face_Normal(1,3,1,1), Face_Normal(2,3,1,1)) + Ewb*AM(U(:,1,1), Ugbw, Face_Normal(1,3,1,1), Face_Normal(2,3,1,1)))*S(3,1,1) ...                       % Bottom (Wall BC)
                     + AP(U(:,1,1), U(:,1,2), Face_Normal(1,4,1,1), Face_Normal(2,4,1,1))*S(4,1,1);                                                         % Top 
    
    f =  -( AP(U(:,1,1), Ugi, Face_Normal(1,1,1,1), Face_Normal(2,1,1,1))*S(1,1,1)*U(:,1,1) + AM(U(:,1,1), Ugi, Face_Normal(1,1,1,1), Face_Normal(2,1,1,1))*S(1,1,1)*Ugi...                           % Left (Inlet BC at time n)
            + AP(U(:,1,1), U(:,2,1), Face_Normal(1,2,1,1), Face_Normal(2,2,1,1))*S(2,1,1)*U(:,1,1) + AM(U(:,1,1), U(:,2,1), Face_Normal(1,2,1,1), Face_Normal(2,2,1,1))*S(2,1,1)*U(:,2,1)...            % Right (time n)
            + AP(U(:,1,1), Ugbw, Face_Normal(1,3,1,1), Face_Normal(2,3,1,1))*S(3,1,1)*U(:,1,1) + AM(U(:,1,1), Ugbw, Face_Normal(1,3,1,1), Face_Normal(2,3,1,1))*S(3,1,1)*Ugbw...                        % Bottom (time n) (Wall BC at time n)
            + AP(U(:,1,1), U(:,1,2), Face_Normal(1,4,1,1), Face_Normal(2,4,1,1))*S(4,1,1)*U(:,1,1) + AM(U(:,1,1), U(:,1,2), Face_Normal(1,4,1,1), Face_Normal(2,4,1,1))*S(4,1,1)*U(:,1,2) )...          % Top      
            - AM(U(:,1,1), U(:,2,1), Face_Normal(1,2,1,1), Face_Normal(2,2,1,1))*S(2,1,1)*dU_old(:,2,1);                                                    % Right delta U at time n
    
    alpha = a - b*g(:,:,2);                 
    v(:,1,1) = alpha\(f-b*v(:,1,2));
    
    dU_new(:,1) = v(:,1,1);
    
    for j = 2:Ny
        dU_new(:,j) = v(:,1,j) - g(:,:,j)*dU_new(:,j-1);
    end

end

%% Function that relaxes the inner lines ( 1 < i < Nx)
function[dU_new] = Inner_Tridiagonal_Solver(dU_old, U, Ny, S, V, Face_Normal, dt, i)
    
    g = zeros(4,4,Ny); 
    v = zeros(4,1,Ny);

    % Top Boundary

    [Ugtw, Ewt] = Wallghost(U(:,i,Ny), Face_Normal(1,4,i,Ny), Face_Normal(2,4,i,Ny));

    a = V(i,Ny)/dt*eye(4) ...
                    + AP(U(:,i,Ny), U(:,i-1,Ny), Face_Normal(1,1,i,Ny), Face_Normal(2,1,i,Ny))*S(1,i,Ny) ...                                                    % Left
                    + AP(U(:,i,Ny), U(:,i+1,Ny), Face_Normal(1,2,i,Ny), Face_Normal(2,2,i,Ny))*S(2,i,Ny) ...                                                     % Right
                    + AP(U(:,i,Ny), U(:,i,Ny-1), Face_Normal(1,3,i,Ny), Face_Normal(2,3,i,Ny))*S(3,i,Ny) ...                                                    % Bottom
                    + (AP(U(:,i,Ny), Ugtw, Face_Normal(1,4,i,Ny), Face_Normal(2,4,i,Ny)) + Ewt*AM(U(:,i,Ny), Ugtw, Face_Normal(1,4,i,Ny), Face_Normal(2,4,i,Ny)))*S(4,i,Ny);                             % Top (Wall BC)
    
    c = AM(U(:,i,Ny), U(:,i,Ny-1), Face_Normal(1,3,i,Ny), Face_Normal(2,3,i,Ny))*S(3,i,Ny); 
    
    f =  -( AP(U(:,i,Ny), U(:,i-1,Ny), Face_Normal(1,1,i,Ny), Face_Normal(2,1,i,Ny))*S(1,i,Ny)*U(:,i,Ny) + AM(U(:,i,Ny), U(:,i-1,Ny), Face_Normal(1,1,i,Ny), Face_Normal(2,1,i,Ny))*S(1,i,Ny)*U(:,i-1,Ny)...        % Left (time n)
            + AP(U(:,i,Ny), U(:,i+1,Ny), Face_Normal(1,2,i,Ny), Face_Normal(2,2,i,Ny))*S(2,i,Ny)*U(:,i,Ny) + AM(U(:,i,Ny), U(:,i+1,Ny), Face_Normal(1,2,i,Ny), Face_Normal(2,2,i,Ny))*S(2,i,Ny)*U(:,i+1,Ny)...        % Right (time n)
            + AP(U(:,i,Ny), U(:,i,Ny-1), Face_Normal(1,3,i,Ny), Face_Normal(2,3,i,Ny))*S(3,i,Ny)*U(:,i,Ny) + AM(U(:,i,Ny), U(:,i,Ny-1), Face_Normal(1,3,i,Ny), Face_Normal(2,3,i,Ny))*S(3,i,Ny)*U(:,i,Ny-1)...      % Bottom (time n)
            + AP(U(:,i,Ny), Ugtw, Face_Normal(1,4,i,Ny), Face_Normal(2,4,i,Ny))*S(4,i,Ny)*U(:,i,Ny) + AM(U(:,i,Ny), Ugtw, Face_Normal(1,4,i,Ny), Face_Normal(2,4,i,Ny))*S(4,i,Ny)*Ugtw )...                           % Top (Wall BC at time n)              
            - AM(U(:,i,Ny), U(:,i+1,Ny), Face_Normal(1,2,i,Ny), Face_Normal(2,2,i,Ny))*S(2,i,Ny)*dU_old(:,i+1,Ny) - AM(U(:,i,Ny), U(:,i-1,Ny), Face_Normal(1,1,i,Ny), Face_Normal(2,1,i,Ny))*S(1,i,Ny)*dU_old(:,i-1,Ny); % Left (time n) and right (time n-1) delta U 
    
    
    alpha = a;
    g(:,:,Ny) = alpha\c;
    v(:,1,Ny) = alpha\f;

    % Middle

    for j = (Ny-1):-1:2

        b = AM(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*S(4,i,j); 

        a = V(i,j)/dt*eye(4) ...
                    + AP(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*S(1,i,j) ...                                                    % Left
                    + AP(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*S(2,i,j) ...                                                     % Right
                    + AP(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*S(3,i,j) ...                                                    % Bottom
                    + AP(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*S(4,i,j);                                                        % Top 
    
        c = AM(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*S(3,i,j); 
    
        f =  -( AP(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*S(1,i,j)*U(:,i,j) + AM(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*S(1,i,j)*U(:,i-1,j)...        % Left (time n)
                + AP(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*S(2,i,j)*U(:,i,j) + AM(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*S(2,i,j)*U(:,i+1,j)...        % Right (time n)
                + AP(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*S(3,i,j)*U(:,i,j) + AM(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*S(3,i,j)*U(:,i,j-1)...      % Bottom (time n)
                + AP(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*S(4,i,j)*U(:,i,j) + AM(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*S(4,i,j)*U(:,i,j+1) )...      % Top (time n)               
                - AM(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*S(2,i,j)*dU_old(:,i+1,j) - AM(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*S(1,i,j)*dU_old(:,i-1,j);   % Left (time n) and right (time n-1) delta U 
    
        alpha = a - b*g(:,:,j+1) ;
        g(:,:,j) = alpha\c;
        v(:,1,j) = alpha\(f-b*v(:,1,j+1));
    end
    
    % Bottom Boundary

    [Ugbw, Ewb] = Wallghost(U(:,i,1), Face_Normal(1,3,i,1), Face_Normal(2,3,i,1));
    
    b = AM(U(:,i,1), U(:,i,2), Face_Normal(1,4,i,1), Face_Normal(2,4,i,1))*S(4,i,1); 
    
    a = V(i,1)/dt*eye(4) ... 
                + AP(U(:,i,1), U(:,i-1,1), Face_Normal(1,1,i,1), Face_Normal(2,1,i,1))*S(1,i,1) ...                                                    % Left
                + AP(U(:,i,1), U(:,i+1,1), Face_Normal(1,2,i,1), Face_Normal(2,2,i,1))*S(2,i,1) ...                                                     % Right
                + (AP(U(:,i,1), Ugbw, Face_Normal(1,3,i,1), Face_Normal(2,3,i,1)) + Ewb*AM(U(:,i,1), Ugbw, Face_Normal(1,3,i,1), Face_Normal(2,3,i,1)))*S(3,i,1) ...                        % Bottom (Wall BC)
                + AP(U(:,i,1), U(:,i,2), Face_Normal(1,4,i,1), Face_Normal(2,4,i,1))*S(3,i,1);                                                          % Top 
    
    f = - AP(U(:,i,1), U(:,i-1,1), Face_Normal(1,1,i,1), Face_Normal(2,1,i,1))*S(1,i,1)*U(:,i,1) - AM(U(:,i,1), U(:,i-1,1), Face_Normal(1,1,i,1), Face_Normal(2,1,i,1))*S(1,i,1)*U(:,i-1,1)...        % Left (time n)
        - AP(U(:,i,1), U(:,i+1,1), Face_Normal(1,2,i,1), Face_Normal(2,2,i,1))*S(2,i,1)*U(:,i,1) - AM(U(:,i,1), U(:,i+1,1), Face_Normal(1,2,i,1), Face_Normal(2,2,i,1))*S(2,i,1)*U(:,i+1,1)...        % Right (time n)
        - AP(U(:,i,1), Ugbw, Face_Normal(1,3,i,1), Face_Normal(2,3,i,1))*S(3,i,1)*U(:,i,1) - AM(U(:,i,1), Ugbw, Face_Normal(1,3,i,1), Face_Normal(2,3,i,1))*S(3,i,1)*Ugbw...                        % Bottom (Wall BC at time n)
        - AP(U(:,i,1), U(:,i,2), Face_Normal(1,4,i,1), Face_Normal(2,4,i,1))*S(4,i,1)*U(:,i,1) - AM(U(:,i,1), U(:,i,2), Face_Normal(1,4,i,1), Face_Normal(2,4,i,1))*S(4,i,1)*U(:,i,2)...            % Top (time n)            
        - AM(U(:,i,1), U(:,i+1,1), Face_Normal(1,2,i,1), Face_Normal(2,2,i,1))*S(2,i,1)*dU_old(:,i+1,1) - AM(U(:,i,1), U(:,i-1,1), Face_Normal(1,1,i,1), Face_Normal(2,1,i,1))*S(1,i,1)*dU_old(:,i-1,1);   % Left (time n) and right (time n-1) delta U 
    
    
    alpha = a - b*g(:,:,2);                
    v(:,1,1) = alpha\(f-b*v(:,1,2));

    dU_new(:,1) = v(:,1,1);
    
    for j = 2:Ny
        dU_new(:,j) = v(:,1,j) - g(:,:,j)*dU_new(:,j-1);
    end

end



%% Function that relaxes the right-most line
function[dU] = Right_Boundary_Tridiagonal_Solver(dU_old, U, Ny, Nx, S, V, Face_Normal, dt)
   
    g = zeros(4,4,Ny); 
    v = zeros(4,1,Ny);
    Eo = eye(4);    % Outlet E matrix is just the identity matrix    
    
    % Top
    
    [Ugtw, Ewt] = Wallghost(U(:,Nx,Ny), Face_Normal(1,4,Nx,Ny), Face_Normal(2,4,Nx,Ny));
    Ugo = U(:,Nx,Ny);

    a = V(Nx,Ny)/dt*eye(4)...
                    + AP(U(:,Nx,Ny), U(:,Nx-1,Ny), Face_Normal(1,1,Nx,Ny), Face_Normal(2,1,Nx,Ny))*S(1,Nx,Ny) ...                                                    % Left
                    + (AP(U(:,Nx,Ny), Ugo, Face_Normal(1,2,Nx,Ny), Face_Normal(2,2,Nx,Ny)) + Eo*AM(U(:,Nx,Ny), Ugo, Face_Normal(1,2,Nx,Ny), Face_Normal(2,2,Nx,Ny)))*S(2,Nx,Ny)...                               % Right (Outflow BC)
                    + AP(U(:,Nx,Ny), U(:,Nx,Ny-1), Face_Normal(1,3,Nx,Ny), Face_Normal(2,3,Nx,Ny))*S(3,Nx,Ny) ...                                                    % Bottom
                    + (AP(U(:,Nx,Ny), Ugtw, Face_Normal(1,4,Nx,Ny), Face_Normal(2,4,Nx,Ny)) + Ewt*AM(U(:,Nx,Ny), Ugtw, Face_Normal(1,4,Nx,Ny), Face_Normal(2,4,Nx,Ny)))*S(4,Nx,Ny);                             % Top (Wall BC)
    
    c = AM(U(:,Nx,Ny), U(:,Nx,Ny-1), Face_Normal(1,3,Nx,Ny), Face_Normal(2,3,Nx,Ny))*S(3,Nx,Ny); % Correct
    
    f =  -( AP(U(:,Nx,Ny), U(:,Nx-1,Ny), Face_Normal(1,1,Nx,Ny), Face_Normal(2,1,Nx,Ny))*S(1,Nx,Ny)*U(:,Nx,Ny) + AM(U(:,Nx,Ny), U(:,Nx-1,Ny), Face_Normal(1,1,Nx,Ny), Face_Normal(2,1,Nx,Ny))*S(1,Nx,Ny)*U(:,Nx-1,Ny)...        % Left (time n)                       
            + AP(U(:,Nx,Ny), Ugo, Face_Normal(1,1,Nx,Ny), Face_Normal(2,1,Nx,Ny))*S(2,Nx,Ny)*U(:,Nx,Ny) + AM(U(:,Nx,Ny), Ugo, Face_Normal(1,2,Nx,Ny), Face_Normal(2,2,Nx,Ny))*S(2,Nx,Ny)*Ugo...                             % Right (Outlet BC at time n)
            + AP(U(:,Nx,Ny), U(:,Nx,Ny-1), Face_Normal(1,3,Nx,Ny), Face_Normal(2,3,Nx,Ny))*S(3,Nx,Ny)*U(:,Nx,Ny) + AM(U(:,Nx,Ny), U(:,Nx,Ny-1), Face_Normal(1,3,Nx,Ny), Face_Normal(2,3,Nx,Ny))*S(3,Nx,Ny)*U(:,Nx,Ny-1)...      % Bottom (time n)
            + AP(U(:,Nx,Ny), Ugtw, Face_Normal(1,4,Nx,Ny), Face_Normal(2,4,Nx,Ny))*S(4,Nx,Ny)*U(:,Nx,Ny) + AM(U(:,Nx,Ny), Ugtw, Face_Normal(1,4,Nx,Ny), Face_Normal(2,4,Nx,Ny))*S(4,Nx,Ny)*Ugtw )...                       % Top (Wall BC at time n)               
            - AM(U(:,Nx,Ny), U(:,Nx-1,Ny), Face_Normal(1,1,Nx,Ny), Face_Normal(2,1,Nx,Ny))*S(1,Nx,Ny)*dU_old(:,Nx-1,Ny);                                                      % Left delta U at time n
    
    
    alpha = a;
    g(:,:,Ny) = alpha\c;
    v(:,1,Ny) = alpha\f;
    
    % Middle

    for j = (Ny-1):-1:2    
        
        Ugo = U(:,Nx,j);
        b = AM( U(:,Nx,j), U(:,Nx,j+1), Face_Normal(1,4,Nx,j), Face_Normal(2,4,Nx,j))*S(4,Nx,j); 
    
        a = V(Nx,j)/dt*eye(4) ...
                + AP(U(:,Nx,j), U(:,Nx-1,j), Face_Normal(1,1,Nx,j), Face_Normal(2,1,Nx,j))*S(1,Nx,j) ...                                                    % Left
                + (AP(U(:,Nx,j), Ugo, Face_Normal(1,2,Nx,j), Face_Normal(2,2,Nx,j)) + Eo*AM(U(:,Nx,j), Ugo, Face_Normal(1,2,Nx,j), Face_Normal(2,2,Nx,j)))*S(2,Nx,j)...                                % Right (Outflow BC)
                + AP(U(:,Nx,j), U(:,Nx,j-1), Face_Normal(1,3,Nx,j), Face_Normal(2,3,Nx,j))*S(3,Nx,j) ...                                                    % Bottom
                + AP(U(:,Nx,j), U(:,Nx,j+1), Face_Normal(1,4,Nx,j), Face_Normal(2,4,Nx,j))*S(4,Nx,j);                                                        % Top
    
        c = AM(U(:,Nx,j), U(:,Nx,j-1), Face_Normal(1,3,Nx,j), Face_Normal(2,3,Nx,j))*S(3,Nx,j); 
    
        f =  -( AP(U(:,Nx,j), U(:,Nx-1,j), Face_Normal(1,1,Nx,j), Face_Normal(2,1,Nx,j))*S(1,Nx,j)*U(:,Nx,j) + AM(U(:,Nx,j), U(:,Nx-1,j), Face_Normal(1,1,Nx,j), Face_Normal(2,1,Nx,j))*S(1,Nx,j)*U(:,Nx-1,j)...        % Left (time n)                       
                + AP(U(:,Nx,j), Ugo, Face_Normal(1,2,Nx,j), Face_Normal(2,2,Nx,j))*S(2,Nx,j)*U(:,Nx,j) + AM(U(:,Nx,j), Ugo, Face_Normal(1,2,Nx,j), Face_Normal(2,2,Nx,j))*S(2,Nx,j)*Ugo...                             % Right (Outlet BC at time n)
                + AP(U(:,Nx,j), U(:,Nx,j-1), Face_Normal(1,3,Nx,j), Face_Normal(2,3,Nx,j))*S(3,Nx,j)*U(:,Nx,j) + AM(U(:,Nx,j), U(:,Nx,j-1), Face_Normal(1,3,Nx,j), Face_Normal(2,3,Nx,j))*S(3,Nx,j)*U(:,Nx,j-1)...      % Bottom (time n)
                + AP(U(:,Nx,j), U(:,Nx,j+1), Face_Normal(1,4,Nx,j), Face_Normal(2,4,Nx,j))*S(4,Nx,j)*U(:,Nx,j) + AM(U(:,Nx,j), U(:,Nx,j+1), Face_Normal(1,4,Nx,j), Face_Normal(2,4,Nx,j))*S(4,Nx,j)*U(:,Nx,j+1) )...      % Top (time n)               
                - AM(U(:,Nx,j), U(:,Nx-1,j), Face_Normal(1,1,Nx,j), Face_Normal(2,1,Nx,j))*S(1,Nx,j)*dU_old(:,Nx-1,j);                                                     % Left delta U at time n
    
        alpha = a - b*g(:,:,j+1);
        g(:,:,j) = alpha\c;
        v(:,1,j) = alpha\(f-b*v(:,1,j+1));
    
    end
    
    % Bottom

    [Ugbw, Ewb] = Wallghost(U(:,Nx,1), Face_Normal(1,3,Nx,1), Face_Normal(2,3,Nx,1));
    Ugo = U(:,Nx,1);
    
    b = AM(U(:,Nx,1), U(:,Nx,2), Face_Normal(1,4,Nx,1), Face_Normal(2,4,Nx,1))*S(4,Nx,1); 
    
    a = V(Nx,1)/dt*eye(4) ...
                    + AP(U(:,Nx,1), U(:,Nx-1,1), Face_Normal(1,1,Nx,1), Face_Normal(2,1,Nx,1))*S(1,Nx,1) ...                                                    % Left
                    + (AP(U(:,Nx,1), Ugo, Face_Normal(1,2,Nx,1), Face_Normal(2,2,Nx,1)) + Eo*AM(U(:,Nx,1), Ugo, Face_Normal(1,2,Nx,1), Face_Normal(2,2,Nx,1)))*S(2,Nx,1)...                              % Right (Outflow BC)
                    + (AP(U(:,Nx,1), Ugbw, Face_Normal(1,3,Nx,1), Face_Normal(2,3,Nx,1)) + Ewb*AM(U(:,Nx,1), Ugbw, Face_Normal(1,3,Nx,1), Face_Normal(2,3,Nx,1)))*S(3,Nx,1) ...                        % Bottom (Wall BC) 
                    + AP(U(:,Nx,1), U(:,Nx,2), Face_Normal(1,4,Nx,1), Face_Normal(2,4,Nx,1))*S(4,Nx,1);                                                          % Top  
    
    f =  -( AP(U(:,Nx,1), U(:,Nx-1,1), Face_Normal(1,1,Nx,1), Face_Normal(2,1,Nx,1))*S(1,Nx,1)*U(:,Nx,1) + AM(U(:,Nx,1), U(:,Nx-1,1), Face_Normal(1,1,Nx,1), Face_Normal(2,1,Nx,1))*S(1,Nx,1)*U(:,Nx-1,1)...        % Left (time n)                       
            + AP(U(:,Nx,1), Ugo, Face_Normal(1,2,Nx,1), Face_Normal(2,2,Nx,1))*S(2,Nx,1)*U(:,Nx,1) + AM(U(:,Nx,1), Ugo, Face_Normal(1,2,Nx,1), Face_Normal(2,2,Nx,1))*S(2,Nx,1)*Ugo...                             % Right (Outlet BC at time n)
            + AP(U(:,Nx,1), Ugbw, Face_Normal(1,3,Nx,1), Face_Normal(2,3,Nx,1))*S(3,Nx,1)*U(:,Nx,1) + AM(U(:,Nx,1), Ugbw, Face_Normal(1,3,Nx,1), Face_Normal(2,3,Nx,1))*S(3,Nx,1)*Ugbw...                       % Bottom (Wall BC at time n)
            + AP(U(:,Nx,1), U(:,Nx,2), Face_Normal(1,4,Nx,1), Face_Normal(2,4,Nx,1))*S(4,Nx,1)*U(:,Nx,1) + AM(U(:,Nx,1), U(:,Nx,2), Face_Normal(1,4,Nx,1), Face_Normal(2,4,Nx,1))*S(4,Nx,1)*U(:,Nx,2) )...            % Top (time n)              
            - AM(U(:,Nx,1), U(:,Nx-1,1), Face_Normal(1,1,Nx,1), Face_Normal(2,1,Nx,1))*S(1,Nx,1)*dU_old(:,Nx-1,1);                                                     % Left delta U at time n
    
    
    alpha = a - b*g(:,:,2);                
    v(:,1,1) = alpha\(f-b*v(:,1,2));
   
    dU(:,1) = v(:,1,1);
    
    for j = 2:Ny
        dU(:,j) = v(:,1,j) - g(:,:,j)*dU(:,j-1);
    end

end


%% Function for positive eigenvalue flux Jacobian

function [Aplus] = AP(Ui, Uii, nx, ny)

gamma = 1.4;

Wi = consToPrim(Ui);
Wii = consToPrim(Uii);

dp = abs(Wi(4) - Wii(4))/min(Wi(4),Wii(4));
g = 5.72;
weight = 1 - 0.5*(1/((g*dp)^2 + 1));

Up = weight*Ui + (1-weight)*Uii;
Wp = consToPrim(Up);

rho = Wp(1);
u = Wp(2);
v = Wp(3);
P = Wp(4);

a = sqrt(gamma*P/rho);

uprime = u*nx + v*ny;
l = [1/2*(uprime - a + abs(uprime - a)), 1/2*(uprime + abs(uprime)), 1/2*(uprime + abs(uprime)), 1/2*(uprime + a + abs(uprime + a))];

M = [l(2), rho*nx/(2*a)*(-l(1) + l(4)), rho*ny/(2*a)*(-l(1) + l(4)), 1/(2*a^2)*(l(1) - 2*l(2) + l(4));
    0, 1/2*(l(1)*nx^2 + 2*l(3)*ny^2 + l(4)*nx^2), nx*ny/2*(l(1) - 2*l(3) + l(4)), nx/(2*a*rho)*(-l(1) + l(4));
    0, nx*ny/2*(l(1) - 2*l(3) + l(4)), 1/2*(l(1)*ny^2 + 2*l(3)*nx^2 + l(4)*ny^2), ny/(2*a*rho)*(-l(1) + l(4));
    0, a*rho*nx/2*(-l(1) + l(4)), a*rho*ny/2*(-l(1) + l(4)), 1/2*(l(1) + l(4))];


dvdu = [1, 0, 0, 0;
    -u/rho, 1/rho, 0, 0;
    -v/rho 0, 1/rho, 0;
    1/2*(gamma-1)*(u^2 + v^2), -(gamma-1)*u, -(gamma-1)*v, gamma - 1];

dudv = [1, 0, 0, 0;
    u, rho, 0, 0;
    v, 0, rho, 0;
    1/2*(u^2+v^2), rho*u, rho*v, 1/(gamma-1)];


Aplus = dudv*M*dvdu;

end

%% Function for negative eigenvalue flux Jacobian

function [Aminus] = AM(Ui, Uii, nx, ny)
gamma = 1.4;

Wi = consToPrim(Ui);
Wii = consToPrim(Uii);

dp = abs(Wi(4) - Wii(4))/min(Wi(4),Wii(4));
g = 5.72;
weight = 1 - 0.5*(1/((g*dp)^2 + 1));

Um = (1-weight)*Ui + weight*Uii;
Wm = consToPrim(Um);

rho = Wm(1);
u = Wm(2);
v = Wm(3);
P = Wm(4);

a = sqrt(gamma*P/rho);

uprime = u*nx + v*ny;
l = [1/2*(uprime - a - abs(uprime - a)), 1/2*(uprime - abs(uprime)), 1/2*(uprime - abs(uprime)), 1/2*(uprime + a - abs(uprime + a))];

M = [l(2), rho*nx/(2*a)*(-l(1) + l(4)), rho*ny/(2*a)*(-l(1) + l(4)), 1/(2*a^2)*(l(1) - 2*l(2) + l(4));
    0, 1/2*(l(1)*nx^2 + 2*l(3)*ny^2 + l(4)*nx^2), nx*ny/2*(l(1) - 2*l(3) + l(4)), nx/(2*a*rho)*(-l(1) + l(4));
    0, nx*ny/2*(l(1) - 2*l(3) + l(4)), 1/2*(l(1)*ny^2 + 2*l(3)*nx^2 + l(4)*ny^2), ny/(2*a*rho)*(-l(1) + l(4));
    0, a*rho*nx/2*(-l(1) + l(4)), a*rho*ny/2*(-l(1) + l(4)), 1/2*(l(1) + l(4))];

dvdu = [1, 0, 0, 0;
    -u/rho, 1/rho, 0, 0;
    -v/rho 0, 1/rho, 0;
    1/2*(gamma-1)*(u^2 + v^2), -(gamma-1)*u, -(gamma-1)*v, gamma - 1];

dudv = [1, 0, 0, 0;
    u, rho, 0, 0;
    v, 0, rho, 0;
    1/2*(u^2+v^2), rho*u, rho*v, 1/(gamma-1)];

Aminus = dudv*M*dvdu;

end

%% Function that converts primitive variables to conserved variables

function [U] = primToCons(W)

    gamma = 1.4;

    % Extract primitive variables from W
    rho = W(1);  % Density
    u   = W(2);  % Velocity in x-direction
    v   = W(3);  % Velocity in y-direction
    p   = W(4);  % Pressure

    % Calculate conserved variables
    U = zeros(4, 1);
    U(1) = rho;               % rho
    U(2) = rho * u;          % rho * u
    U(3) = rho * v;          % rho * v
    U(4) = p / (gamma - 1) + 0.5 * rho * (u^2 + v^2);  % Total energy

end

%% Function that converts conserved variables to primitive variables

function [W] = consToPrim(U)
           
    gamma = 1.4;    

    % Extract conserved variables from U
    rho = U(1);          % Density
    rhou = U(2);         % rho * u
    rhov = U(3);         % rho * v
    E = U(4);            % Total energy

    % Calculate primitive variables
    W = zeros(4,1);
    W(1) = rho;                  % Density is directly available
    W(2) = rhou / rho;          % Velocity in x-direction: u = (rho * u) / rho
    W(3) = rhov / rho;          % Velocity in y-direction: v = (rho * v) / rho

    % Total energy equation to calculate pressure
    kinetic_energy = 0.5 * rho * (W(2)^2 + W(3)^2);
    W(4) = (gamma - 1) * (E - kinetic_energy);  % Pressure: p = (gamma - 1) * (E - 0.5*rho*(u^2 + v^2))

end
%% Function that calculates U for ghost cells at a wall
function [Ugw,E] = Wallghost(U, nx, ny)
    
    Ghost_primitive = zeros(4,1);
    Inside_Primitives = consToPrim(U);        
    
    u = Inside_Primitives(2);
    v = Inside_Primitives(3);
    
    Ghost_primitive(1) = Inside_Primitives(1);
    Ghost_primitive(2) = u - 2*(u*nx + v*ny)*nx;
    Ghost_primitive(3) = v - 2*(u*nx + v*ny)*ny;
    Ghost_primitive(4) = Inside_Primitives(4);

    Ugw = primToCons(Ghost_primitive);

    E = [1, 0, 0, 0;
         0, (1 - 2*nx^2), -2*nx*ny, 0;
         0, -2*nx*ny,  (1-2*ny^2), 0;
         0, 0, 0, 1];    

end

%% Function that calculates residual for GSLR
function [residual] = calculateLocalResidual(dU, dU_old, U, V, dt, Ny, S, Face_Normal, i)
    
    lineresidual = 0;


    for j = 2:Ny-1

        b = AM(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*S(4,i,j); 

        a = V(i,j)/dt*eye(4) ...
                    + AP(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*S(1,i,j) ...                                                    % Left
                    + AP(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*S(2,i,j) ...                                                     % Right
                    + AP(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*S(3,i,j) ...                                                    % Bottom
                    + AP(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*S(4,i,j);                                                        % Top 
    
        c = AM(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*S(3,i,j); 
    
        f =  -( AP(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*S(1,i,j)*U(:,i,j) + AM(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*S(1,i,j)*U(:,i-1,j)...        % Left (time n)
                + AP(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*S(2,i,j)*U(:,i,j) + AM(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*S(2,i,j)*U(:,i+1,j)...        % Right (time n)
                + AP(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*S(3,i,j)*U(:,i,j) + AM(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*S(3,i,j)*U(:,i,j-1)...      % Bottom (time n)
                + AP(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*S(4,i,j)*U(:,i,j) + AM(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*S(4,i,j)*U(:,i,j+1) )...      % Top (time n)               
                - AM(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*S(2,i,j)*dU_old(:,i+1,j) - AM(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*S(1,i,j)*dU_old(:,i-1,j);   % Left (time n) and right (time n-1) delta U 
    
                res =  b*dU(:,i,j+1) + a*dU(:,i,j) + c*dU(:,i,j-1) - f;      
                lineresidual = lineresidual + res(1)^2;

    end

    residual = sqrt(lineresidual);

end

%% Function that calculated outer residual
function [total_residual] = Calculate_Residual(U, V, S, Face_Normal, Nx, Ny)

    residual = 0;
    
    for i=2:Nx-1
        for j = 2:Ny-1
    
            residual =  residual + ((AP(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*S(1,i,j)*U(:,i,j) + AM(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*S(1,i,j)*U(:,i-1,j)...        % Left (time n)
                + AP(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*S(2,i,j)*U(:,i,j) + AM(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*S(2,i,j)*U(:,i+1,j)...        % Right (time n)
                + AP(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*S(3,i,j)*U(:,i,j) + AM(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*S(3,i,j)*U(:,i,j-1)...      % Bottom (time n)
                + AP(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*S(4,i,j)*U(:,i,j) + AM(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*S(4,i,j)*U(:,i,j+1))/V(i,j)).^2;      % Top (time n)               
               
        end
    end
   
    res = residual(1);
    total_residual = sqrt(res);

end

%% Function that calculates dt
function [dt] = Calculate_dt(U, S, Nx, Ny,CFL)

    gamma = 1.4;
    int_dt = zeros(Nx,Ny);
    
    for i = 1:Nx
        for j = 1:Ny

            V = consToPrim(U(:,i,j));
            dx = min(min(S(1,i,j), S(2,i,j)));
            dy = min(min(S(3,i,j), S(4,i,j)));
            c = sqrt(gamma*V(4)/V(1));

            int_dt(i,j) =  CFL/( abs(V(2)/dx) + abs(V(3)/dy) + c*sqrt(1/dx^2 + 1/dy^2) );

        end
    end

    dt = min(min(int_dt));
    
end

