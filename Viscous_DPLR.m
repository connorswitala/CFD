clear all;
close all;
clc;

%% Grid geometry and creation

% run('CylinderGridGenerator.m');
run('Flat_Plate_Generator.m');

%% Inflow conditions

% M = 5;
% T = 300;
% gamma = 1.4;
% R = 287;
% a = sqrt(gamma*R*T);
% P = 10000;
% rho = P/(R*T);
% u = M*a;
% v = 0;

M = 1.5;
T = 300;
gamma = 1.4;
R = 287;
a = sqrt(gamma*R*T);
rho = 0.01;
P = rho*R*T;
u = M*a;
v = 0;


%% Initialize U

W_inlet = [rho; u; v; P];
U_inlet = primToCons(W_inlet);
U_old = zeros(4,Nx,Ny);

xc = squeeze(center(1,:,:));
yc = squeeze(center(2,:,:));


% for j = 1:Ny
%     for i = 1:Nx
%         U_old(:,i,j) = U_inlet;
%     end
% end


%% Start physics solving
% 
% dU = zeros(4,Nx,Ny); % Initial guess for dU
% t_old = 0;
% outer_residual = 10;
% n = 1;
% U = U_old;

load('100x100_viscous_plate_simulation_data.mat', 'U', 'dU', 't', 'n', 'outer_residual');
U_old = U;
t_old = t(n-1);

dt = Calculate_dt(U, S, Nx, Ny);

while outer_residual >= 10e-9
    
    dU_old = dU;
    U_old = U;
    dU = SolveOneTimestep(dU_old, U_old, U_inlet, Nx, Ny, S, V, Face_Normal, dt, xc, yc);
    U = U_old + dU;    
    t(n) = t_old + dt;
    t_old = t(n);
    outer_residual(n) = Calculate_Residual(U, V, S, Face_Normal, Nx, Ny, xc, yc);

    if abs(mod(n, 20)) == 0
        run('Postprocessingplate.m');
        % save('100x100_viscous_plate_simulation_data.mat', 'U', 'dU', 't', 'n', 'outer_residual')
        look = sprintf('Outer Residual = %e,  dt = %e, Iteration number = %d.\n', outer_residual(n), dt, n);
        disp(look);
    end
    

    dt = Calculate_dt(U, S, Nx, Ny);
    n = n+1;
    

end



%% Function that goes through each lines and calls the propper GSLR for it. Repeats until convergence.

function [dU] = SolveOneTimestep(dU_old, U, U_inlet, Nx, Ny, S, V, Face_Normal, dt, xc, yc)
    
    global_inner_residual = 1;
    dU_new = zeros(4,Nx,Ny);
  
    while global_inner_residual >= 10e-11

        global_inner_residual = 0;

        for i = 1:Nx
            
            if i == 1

                dU_new(:,1,:) = Left_Boundary_Tridiagonal_Solver(dU_old, U, U_inlet, Ny, S, V, Face_Normal, dt, xc, yc);  

            elseif (i ~= 1 && i ~= Nx)

                dU_new(:,i,:) = Inner_Tridiagonal_Solver(dU_old, U, U_inlet, Ny, S, V, Face_Normal, dt, i, xc, yc);
                global_inner_residual = global_inner_residual + (calculateLocalResidual(dU_new, dU_old, U, V, dt, Ny, S, Face_Normal, i, xc, yc))^2;

            elseif i == Nx

                dU_new(:,Nx,:) = Right_Boundary_Tridiagonal_Solver(dU_old, U, U_inlet, Ny, Nx, S, V, Face_Normal, dt, xc, yc);               

            end

        end

        global_inner_residual = sqrt(global_inner_residual)
        dU_old = dU_new;

    end

    dU = dU_old;

end

%% Function that relaxes the left-most line
function[dU_new] = Left_Boundary_Tridiagonal_Solver(dU_old, U, U_inlet, Ny, S, V, Face_Normal, dt, xc, yc)
    
    g = zeros(4,4,Ny); 
    v = zeros(4,1,Ny);
    dU_new = zeros(4,Ny);
    
    % [UgiT, UgvT, EiT, EvT] = Inlet(U_inlet, Face_Normal(1,4,1,Ny), Face_Normal(2,4,1,Ny));
    % [UgiL, UgvL, EiL, EvL] = Outlet(U(:,1,Ny), Face_Normal(1,1,1,Ny), Face_Normal(2,1,1,Ny)); 

    [UgiT, UgvT, EiT, EvT] = Outlet(U(:,1,Ny), Face_Normal(1,4,1,Ny), Face_Normal(2,4,1,Ny));
    [UgiL, UgvL, EiL, EvL] = Inlet(U_inlet, Face_Normal(1,1,1,Ny), Face_Normal(2,1,1,Ny)); 

    
    % Top

    a =  V(1,Ny)/dt*eye(4) ...
             + ( AP(U(:,1,Ny), UgiL, Face_Normal(1,1,1,Ny), Face_Normal(2,1,1,Ny))...                                                                                                                            % Left (Inlet BC)
             + EiL*AM(U(:,1,Ny), UgiL, Face_Normal(1,1,1,Ny), Face_Normal(2,1,1,Ny))...                                                                                                                          % Left (Inlet BC)
             + ImpV_bx(U(:,1,Ny), UgvL, EvL, xc(1,Ny), xc(2,Ny), Face_Normal(1,1,1,Ny), Face_Normal(2,1,1,Ny)) )*S(1,1,Ny) ...                                                                                   % Left (Inlet BC) 
             + ( AP(U(:,1,Ny), U(:,2,Ny), Face_Normal(1,2,1,Ny), Face_Normal(2,2,1,Ny)) - ImpV_x(U(:,1,Ny), U(:,2,Ny), xc(1,Ny), xc(2,Ny), Face_Normal(1,2,1,Ny), Face_Normal(2,2,1,Ny)) )*S(2,1,Ny) ...        % Right
             + ( AP(U(:,1,Ny), U(:,1,Ny-1), Face_Normal(1,3,1,Ny), Face_Normal(2,3,1,Ny)) + ImpV_y(U(:,1,Ny), U(:,1,Ny-1), yc(1,Ny), yc(1,Ny-1), Face_Normal(1,3,1,Ny), Face_Normal(2,3,1,Ny)) )*S(3,1,Ny) ...  % Bottom
             + ( AP(U(:,1,Ny), UgiT, Face_Normal(1,4,1,Ny), Face_Normal(2,4,1,Ny))...                                                                                                                           % Top (Wall BC)
             + EiT*AM(U(:,1,Ny), UgiT, Face_Normal(1,4,1,Ny), Face_Normal(2,4,1,Ny))...                                                                                                                        % Top (Wall BC)
             + ImpV_by(U(:,1,Ny), UgvT, EvT, yc(1,Ny), yc(1,Ny-1), Face_Normal(1,4,1,Ny), Face_Normal(2,4,1,Ny)) )*S(4,1,Ny);                                                                                 % Top (Wall BC)                
    
    c = ( AM(U(:,1,Ny), U(:,1,Ny-1), Face_Normal(1,3,1,Ny), Face_Normal(2,3,1,Ny)) - ImpV_y(U(:,1,Ny), U(:,1,Ny-1), yc(1,Ny), yc(1,Ny-1), Face_Normal(1,3,1,Ny), Face_Normal(2,3,1,Ny)) )*S(3,1,Ny);            % Bottom delta U
    
    f =  -( ( AP(U(:,1,Ny), UgiL, Face_Normal(1,1,1,Ny), Face_Normal(2,1,1,Ny))*U(:,1,Ny) ...                                                                                                                            % Left (Inlet BC)
            + AM(U(:,1,Ny), UgiL, Face_Normal(1,1,1,Ny), Face_Normal(2,1,1,Ny))*UgiL + ExpV_x(U(:,1,Ny), UgvL, xc(1,Ny), xc(2,Ny), Face_Normal(1,1,1,Ny), Face_Normal(2,1,1,Ny)) )*S(1,1,Ny)...                            % Left (Inlet BC) 
            + ( AP(U(:,1,Ny), U(:,2,Ny), Face_Normal(1,2,1,Ny), Face_Normal(2,2,1,Ny))*U(:,1,Ny)...                                                                                                                     % Right
            + AM(U(:,1,Ny), U(:,2,Ny), Face_Normal(1,2,1,Ny), Face_Normal(2,2,1,Ny))*U(:,2,Ny) + ExpV_x(U(:,1,Ny), U(:,2,Ny), xc(1,Ny), xc(2,Ny), Face_Normal(1,2,1,Ny), Face_Normal(2,2,1,Ny)) )*S(2,1,Ny)...          % Right  
            + ( AP(U(:,1,Ny), U(:,1,Ny-1), Face_Normal(1,3,1,Ny), Face_Normal(2,3,1,Ny))*U(:,1,Ny)...                                                                                                                   % Bottom
            + AM(U(:,1,Ny), U(:,1,Ny-1), Face_Normal(1,3,1,Ny), Face_Normal(2,3,1,Ny))*U(:,1,Ny-1) + ExpV_y(U(:,1,Ny), U(:,1,Ny-1), yc(1,Ny), yc(1,Ny-1), Face_Normal(1,3,1,Ny), Face_Normal(2,3,1,Ny)) )*S(3,1,Ny)...  % Bottom
            + ( AP(U(:,1,Ny), UgiT, Face_Normal(1,4,1,Ny), Face_Normal(2,4,1,Ny))*U(:,1,Ny)...                                                                                                                          % Top (Wall BC)
            + AM(U(:,1,Ny), UgiT, Face_Normal(1,4,1,Ny), Face_Normal(2,4,1,Ny))*UgiT + ExpV_y(U(:,1,Ny), UgvT, yc(1,Ny), yc(1,Ny-1), Face_Normal(1,4,1,Ny), Face_Normal(2,4,1,Ny)) )*S(4,1,Ny) ...                     % Top (Wall BC)
            + ( AM(U(:,1,Ny), U(:,2,Ny), Face_Normal(1,2,1,Ny), Face_Normal(2,2,1,Ny)) + ImpV_x(U(:,1,Ny), U(:,2,Ny), xc(1,Ny), xc(2,Ny), Face_Normal(1,2,1,Ny), Face_Normal(2,2,1,Ny)) )*S(2,1,Ny)*dU_old(:,2,Ny) );   % Right delta U at time n
    
    alpha = a;
    g(:,:,Ny) = alpha\c;
    v(:,1,Ny) = alpha\f;

    % Middle

    for j = (Ny-1):-1:2

        % [UgiL, UgvL, EiL, EvL] = Outlet(U(:,1,j), Face_Normal(1,1,1,j), Face_Normal(2,1,1,j));
        [UgiL, UgvL, EiL, EvL] = Inlet(U_inlet, Face_Normal(1,1,1,j), Face_Normal(2,1,1,j)); 
    

        b = ( AM(U(:,1,j), U(:,1,j+1), Face_Normal(1,4,1,j), Face_Normal(2,4,1,j)) + ImpV_y(U(:,1,j), U(:,1,j+1), yc(1,j), yc(1,j+1), Face_Normal(1,4,1,j), Face_Normal(2,4,1,j)) )*S(4,1,j);               % Top delta U
    
        a =  V(1,j)/dt*eye(4) ...
                 + ( AP(U(:,1,j), UgiL, Face_Normal(1,1,1,j), Face_Normal(2,1,1,j))...                                                                                                                       % Left (Inlet BC)
                 + EiL*AM(U(:,1,j), UgiL, Face_Normal(1,1,1,j), Face_Normal(2,1,1,j)) ...                                                                                                                    % Left (Inlet BC)
                 + ImpV_bx(U(:,1,j), UgvL, EvL, xc(1,j), xc(2,j), Face_Normal(1,1,1,j), Face_Normal(2,1,1,j)) )*S(1,1,j) ...                                                                                 % Left (Inlet BC)
                 + ( AP(U(:,1,j), U(:,2,j), Face_Normal(1,2,1,j), Face_Normal(2,2,1,j)) - ImpV_x(U(:,1,j), U(:,2,j), xc(1,j), xc(2,j), Face_Normal(1,2,1,j), Face_Normal(2,2,1,j)) )*S(2,1,j) ...           % Right
                 + ( AP(U(:,1,j), U(:,1,j-1), Face_Normal(1,3,1,j), Face_Normal(2,3,1,j)) + ImpV_y(U(:,1,j), U(:,1,j-1), yc(1,j), yc(1,j-1), Face_Normal(1,3,1,j), Face_Normal(2,3,1,j)) )*S(3,1,j) ...     % Bottom
                 + ( AP(U(:,1,j), U(:,1,j+1), Face_Normal(1,4,1,j), Face_Normal(2,4,1,j)) - ImpV_y(U(:,1,j), U(:,1,j+1), yc(1,j), yc(1,j+1), Face_Normal(1,4,1,j), Face_Normal(2,4,1,j)) )*S(4,1,j);        % Top 
    
        c = ( AM(U(:,1,j), U(:,1,j-1), Face_Normal(1,3,1,j), Face_Normal(2,3,1,j)) - ImpV_y(U(:,1,j), U(:,1,j-1), yc(1,j), yc(1,j-1), Face_Normal(1,3,1,j), Face_Normal(2,3,1,j)) )*S(3,1,j);               % Bottom delta U
    
        f =  -( ( AP(U(:,1,j), UgiL, Face_Normal(1,1,1,j), Face_Normal(2,1,1,j))*U(:,1,j)...                                                                                                                     % Left (Inlet BC)
                + AM(U(:,1,j), UgiL, Face_Normal(1,1,1,j), Face_Normal(2,1,1,j))*UgiL + ExpV_x(U(:,1,j), UgvL, xc(1,j), xc(2,j), Face_Normal(1,1,1,j), Face_Normal(2,1,1,j)) )*S(1,1,j)...                         % Left (Inlet BC)
                + ( AP(U(:,1,j), U(:,2,j), Face_Normal(1,2,1,j), Face_Normal(2,2,1,j))*U(:,1,j)...                                                                                                              % Right
                + AM(U(:,1,j), U(:,2,j), Face_Normal(1,2,1,j), Face_Normal(2,2,1,j))*U(:,2,j) + ExpV_x(U(:,1,j), U(:,2,j), xc(1,j), xc(2,j), Face_Normal(1,2,1,j), Face_Normal(2,2,1,j)) )*S(2,1,j)...          % Right    
                + ( AP(U(:,1,j), U(:,1,j-1), Face_Normal(1,3,1,j), Face_Normal(2,3,1,j))*U(:,1,j)...                                                                                                            % Bottom
                + AM(U(:,1,j), U(:,1,j-1), Face_Normal(1,3,1,j), Face_Normal(2,3,1,j))*U(:,1,j-1) + ExpV_y(U(:,1,j), U(:,1,j-1), yc(1,j), yc(1,j-1), Face_Normal(1,3,1,j), Face_Normal(2,3,1,j)) )*S(3,1,j)...  % Bottom
                + ( AP(U(:,1,j), U(:,1,j+1), Face_Normal(1,4,1,j), Face_Normal(2,4,1,j))*U(:,1,j)...                                                                                                            % Top
                + AM(U(:,1,j), U(:,1,j+1), Face_Normal(1,4,1,j), Face_Normal(2,4,1,j))*U(:,1,j+1) + ExpV_y(U(:,1,j), U(:,1,j+1), yc(1,j), yc(1,j+1), Face_Normal(1,4,1,j), Face_Normal(2,4,1,j)) )*S(4,1,j) ... % Top     
                + ( AM(U(:,1,j), U(:,2,j), Face_Normal(1,2,1,j), Face_Normal(2,2,1,j)) + ImpV_x(U(:,1,j), U(:,2,j), xc(1,j), xc(2,j), Face_Normal(1,2,1,j), Face_Normal(2,2,1,j)) )*S(2,1,j)*dU_old(:,2,j) );   % Right delta U
        
        alpha = a - b*g(:,:,j+1);
        g(:,:,j) = alpha\c;
        v(:,1,j) = alpha\(f-b*v(:,1,j+1));
 
    end

    % Bottom
    
    [UgiB, UgvB, EiB, EvB] = Viscous_Wall_Ghost(U(:,1,1), Face_Normal(1,3,1,1), Face_Normal(2,3,1,1));
    % [UgiL, UgvL, EiL, EvL] = Outlet(U(:,1,1), Face_Normal(1,1,1,1), Face_Normal(2,1,1,1));
    [UgiL, UgvL, EiL, EvL] = Inlet(U_inlet, Face_Normal(1,1,1,1), Face_Normal(2,1,1,1)); 

    
    b = ( AM(U(:,1,1), U(:,1,2), Face_Normal(1,4,1,1), Face_Normal(2,4,1,1)) + ImpV_y(U(:,1,1), U(:,1,2), yc(1,1), yc(1,2), Face_Normal(1,4,1,1), Face_Normal(2,4,1,1)) )*S(4,1,Ny);        % Top delta U
    
    a =  V(1,1)/dt*eye(4) ... 
             + ( AP(U(:,1,1), UgiL, Face_Normal(1,1,1,1), Face_Normal(2,1,1,1))...                                                                                                               % Left (Inlet BC)
             + EiL*AM(U(:,1,1), UgiL, Face_Normal(1,1,1,1), Face_Normal(2,1,1,1))...                                                                                                             % Left (Inlet BC)
             + ImpV_bx(U(:,1,1), UgvL, EvL, xc(1,1), xc(2,1), Face_Normal(1,1,1,1), Face_Normal(2,1,1,1)) )*S(1,1,1) ...                                                                         % Left (Inlet BC)
             + ( AP(U(:,1,1), U(:,2,1), Face_Normal(1,2,1,1), Face_Normal(2,2,1,1)) - ImpV_x(U(:,1,1), U(:,2,1), xc(1,1), xc(2,1), Face_Normal(1,2,1,1), Face_Normal(2,2,1,1)) )*S(2,1,1) ...   % Right
             + ( AP(U(:,1,1), UgiB, Face_Normal(1,3,1,1), Face_Normal(2,3,1,1))...                                                                                                              % Bottom (Wall BC)
             + EiB*AM(U(:,1,1), UgiB, Face_Normal(1,3,1,1), Face_Normal(2,3,1,1))...                                                                                                           % Bottom (Wall BC)
             + ImpV_by(U(:,1,1), UgvB, EvB, yc(1,1), yc(1,2), Face_Normal(1,3,1,1), Face_Normal(2,3,1,1)) )*S(3,1,1) ...                                                                      % Bottom (Wall BC)
             + ( AP(U(:,1,1), U(:,1,2), Face_Normal(1,4,1,1), Face_Normal(2,4,1,1)) - ImpV_y(U(:,1,1), U(:,1,2), yc(1,1), yc(1,2), Face_Normal(1,4,1,1), Face_Normal(2,4,1,1)) )*S(4,1,1);      % Top 
    
    f =  -( ( AP(U(:,1,1), UgiL, Face_Normal(1,1,1,1), Face_Normal(2,1,1,1))*U(:,1,1)...                                                                                                                     % Left (Inlet BC)
            + AM(U(:,1,1), UgiL, Face_Normal(1,1,1,1), Face_Normal(2,1,1,1))*UgiL + ExpV_x(U(:,1,1), UgvL, xc(1,1), xc(2,1), Face_Normal(1,1,1,1), Face_Normal(2,1,1,1)) )*S(1,1,1)...                         % Left (Inlet BC) 
            + ( AP(U(:,1,1), U(:,2,1), Face_Normal(1,2,1,1), Face_Normal(2,2,1,1))*U(:,1,1)...                                                                                                              % Right
            + AM(U(:,1,1), U(:,2,1), Face_Normal(1,2,1,1), Face_Normal(2,2,1,1))*U(:,2,1) + ExpV_x(U(:,1,1), U(:,2,1), xc(1,1), xc(2,1), Face_Normal(1,2,1,1), Face_Normal(2,2,1,1)) )*S(2,1,1)...          % Right
            + ( AP(U(:,1,1), UgiB, Face_Normal(1,3,1,1), Face_Normal(2,3,1,1))*U(:,1,1)...                                                                                                                  % Bottom (Wall BC)
            + AM(U(:,1,1), UgiB, Face_Normal(1,3,1,1), Face_Normal(2,3,1,1))*UgiB + ExpV_y(U(:,1,1), UgvB, yc(1,1), yc(1,2), Face_Normal(1,3,1,1), Face_Normal(2,3,1,1)) )*S(3,1,1)...                     % Bottom (Wall BC)
            + ( AP(U(:,1,1), U(:,1,2), Face_Normal(1,4,1,1), Face_Normal(2,4,1,1))*U(:,1,1)...                                                                                                              % Top
            + AM(U(:,1,1), U(:,1,2), Face_Normal(1,4,1,1), Face_Normal(2,4,1,1))*U(:,1,2) + ExpV_y(U(:,1,1), U(:,1,2), yc(1,1), yc(1,2), Face_Normal(1,4,1,1), Face_Normal(2,4,1,1)) )*S(4,1,1)...          % Top     
            + ( AM(U(:,1,1), U(:,2,1), Face_Normal(1,2,1,1), Face_Normal(2,2,1,1)) + ImpV_x(U(:,1,1), U(:,2,1), xc(1,1), xc(2,1), Face_Normal(1,2,1,1), Face_Normal(2,2,1,1)) )*S(2,1,1)*dU_old(:,2,1)  );  % Right delta U at time n
    
    alpha = a - b*g(:,:,2);                 
    v(:,1,1) = alpha\(f-b*v(:,1,2));
    
    dU_new(:,1) = v(:,1,1);
    
    for j = 2:Ny
        dU_new(:,j) = v(:,1,j) - g(:,:,j)*dU_new(:,j-1);
    end

end

%% Function that relaxes the inner lines ( 1 < i < Nx)
function[dU_new] = Inner_Tridiagonal_Solver(dU_old, U, U_inlet, Ny, S, V, Face_Normal, dt, i, xc, yc)
    
    g = zeros(4,4,Ny); 
    v = zeros(4,1,Ny);
   
    % Top Boundary

    % [UgiT, UgvT, EiT, EvT] = Inlet(U_inlet, Face_Normal(1,4,i,Ny), Face_Normal(2,4,i,Ny));
    [UgiT, UgvT, EiT, EvT] = Outlet(U(:,i,Ny), Face_Normal(1,4,i,Ny), Face_Normal(2,4,i,Ny));



  
    a = V(i,Ny)/dt*eye(4) ...
            + ( AP(U(:,i,Ny), U(:,i-1,Ny), Face_Normal(1,1,i,Ny), Face_Normal(2,1,i,Ny)) + ImpV_x(U(:,i,Ny), U(:,i-1,Ny), xc(i,Ny), xc(i-1,Ny), Face_Normal(1,1,i,Ny), Face_Normal(2,1,i,Ny)) )*S(1,i,Ny) ...   % Left
            + ( AP(U(:,i,Ny), U(:,i+1,Ny), Face_Normal(1,2,i,Ny), Face_Normal(2,2,i,Ny)) - ImpV_x(U(:,i,Ny), U(:,i+1,Ny), xc(i,Ny), xc(i+1,Ny), Face_Normal(1,2,i,Ny), Face_Normal(2,2,i,Ny)) )*S(2,i,Ny) ...   % Right
            + ( AP(U(:,i,Ny), U(:,i,Ny-1), Face_Normal(1,3,i,Ny), Face_Normal(2,3,i,Ny)) + ImpV_y(U(:,i,Ny), U(:,i,Ny-1), yc(i,Ny), yc(i,Ny-1), Face_Normal(1,3,i,Ny), Face_Normal(2,3,i,Ny)) )*S(3,i,Ny) ...   % Bottom
            + ( AP(U(:,i,Ny), UgiT, Face_Normal(1,4,i,Ny), Face_Normal(2,4,i,Ny))...                                                                                                                            % Top (Wall BC)
            + EiT*AM(U(:,i,Ny), UgiT, Face_Normal(1,4,i,Ny), Face_Normal(2,4,i,Ny))...                                                                                                                         % Top (Wall BC)
            + ImpV_by(U(:,i,Ny), UgvT, EvT, yc(i,Ny), yc(i,Ny-1), Face_Normal(1,4,i,Ny), Face_Normal(2,4,i,Ny)) )*S(4,i,Ny);                                                                                  % Top (Wall BC)
    
    c = ( AM(U(:,i,Ny), U(:,i,Ny-1), Face_Normal(1,3,i,Ny), Face_Normal(2,3,i,Ny)) - ImpV_y(U(:,i,Ny), U(:,i,Ny-1), yc(i,Ny), yc(i,Ny-1), Face_Normal(1,3,i,Ny), Face_Normal(2,3,i,Ny)) )*S(3,i,Ny);            % Bottom delta U
        
    f =  -( ( AP(U(:,i,Ny), U(:,i-1,Ny), Face_Normal(1,1,i,Ny), Face_Normal(2,1,i,Ny))*U(:,i,Ny)...                                                                                                                             % Left 
            + AM(U(:,i,Ny), U(:,i-1,Ny), Face_Normal(1,1,i,Ny), Face_Normal(2,1,i,Ny))*U(:,i-1,Ny) + ExpV_x(U(:,i,Ny), U(:,i-1,Ny), xc(i,Ny), xc(i-1,Ny), Face_Normal(1,1,i,Ny), Face_Normal(2,1,i,Ny)) )*S(1,i,Ny)...          % Left
            + ( AP(U(:,i,Ny), U(:,i+1,Ny), Face_Normal(1,2,i,Ny), Face_Normal(2,2,i,Ny))*U(:,i,Ny)...                                                                                                                           % Right
            + AM(U(:,i,Ny), U(:,i+1,Ny), Face_Normal(1,2,i,Ny), Face_Normal(2,2,i,Ny))*U(:,i+1,Ny) + ExpV_x(U(:,i,Ny), U(:,i+1,Ny), xc(i,Ny), xc(i+1,Ny), Face_Normal(1,2,i,Ny), Face_Normal(2,2,i,Ny)) )*S(2,i,Ny)...          % Right 
            + ( AP(U(:,i,Ny), U(:,i,Ny-1), Face_Normal(1,3,i,Ny), Face_Normal(2,3,i,Ny))*U(:,i,Ny)...                                                                                                                           % Bottom
            + AM(U(:,i,Ny), U(:,i,Ny-1), Face_Normal(1,3,i,Ny), Face_Normal(2,3,i,Ny))*U(:,i,Ny-1) + ExpV_y(U(:,i,Ny), U(:,i,Ny-1), yc(i,Ny), yc(i,Ny-1), Face_Normal(1,3,i,Ny), Face_Normal(2,3,i,Ny)) ) *S(3,i,Ny)...         % Bottom
            + ( AP(U(:,i,Ny), UgiT, Face_Normal(1,4,i,Ny), Face_Normal(2,4,i,Ny))*U(:,i,Ny)...                                                                                                                                  % Top (Wall BC)
            + AM(U(:,i,Ny), UgiT, Face_Normal(1,4,i,Ny), Face_Normal(2,4,i,Ny))*UgiT + ExpV_y(U(:,i,Ny), UgvT, yc(i,Ny), yc(i,Ny-1), Face_Normal(1,4,i,Ny), Face_Normal(2,4,i,Ny)) )*S(4,i,Ny) ...                             % Top (Wall BC)  
            + ( AM(U(:,i,Ny), U(:,i+1,Ny), Face_Normal(1,2,i,Ny), Face_Normal(2,2,i,Ny)) + ImpV_x(U(:,i,Ny), U(:,i+1,Ny), xc(i,Ny), xc(i+1,Ny), Face_Normal(1,2,i,Ny), Face_Normal(2,2,i,Ny)))*S(2,i,Ny)*dU_old(:,i+1,Ny)...    % Right delta U
            + ( AM(U(:,i,Ny), U(:,i-1,Ny), Face_Normal(1,1,i,Ny), Face_Normal(2,1,i,Ny)) - ImpV_x(U(:,i,Ny), U(:,i-1,Ny), xc(i,Ny), xc(i-1,Ny), Face_Normal(1,1,i,Ny), Face_Normal(2,1,i,Ny)))*S(1,i,Ny)*dU_old(:,i-1,Ny) );    % Left delta U   
    
    
    alpha = a;
    g(:,:,Ny) = alpha\c;
    v(:,1,Ny) = alpha\f;

    % Middle

    for j = (Ny-1):-1:2

        b = ( AM(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j)) + ImpV_y(U(:,i,j), U(:,i,j+1), yc(i,j), yc(i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j)) )*S(4,i,j);   % Top delta U

        a = V(i,j)/dt*eye(4) ...
                    + ( AP(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j)) + ImpV_x(U(:,i,j), U(:,i-1,j), xc(i,j), xc(i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j)) )*S(1,i,j) ... % Left
                    + ( AP(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j)) - ImpV_x(U(:,i,j), U(:,i+1,j), xc(i,j), xc(i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j)) )*S(2,i,j) ... % Right
                    + ( AP(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j)) + ImpV_y(U(:,i,j), U(:,i,j-1), yc(i,j), yc(i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j)) )*S(3,i,j) ... % Bottom
                    + ( AP(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j)) - ImpV_y(U(:,i,j), U(:,i,j+1), yc(i,j), yc(i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j)) )*S(4,i,j);    % Top 
    
        c = ( AM(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j)) - ImpV_y(U(:,i,j), U(:,i,j-1), yc(i,j), yc(i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j)) )*S(3,i,j);   % Bottom delta U
    
        f =  -( ( AP(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*U(:,i,j)...                                                                                                                      % Left
                + AM(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*U(:,i-1,j) + ExpV_x(U(:,i,j), U(:,i-1,j), xc(i,j), xc(i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j)) )*S(1,i,j)...          % Left
                + ( AP(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*U(:,i,j)...                                                                                                                    % Right
                + AM(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*U(:,i+1,j) + ExpV_x(U(:,i,j), U(:,i+1,j), xc(i,j), xc(i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j)) )*S(2,i,j)...          % Right
                + ( AP(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*U(:,i,j)...                                                                                                                    % Bottom
                + AM(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*U(:,i,j-1) + ExpV_y(U(:,i,j), U(:,i,j-1), yc(i,j), yc(i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j)) )*S(3,i,j)...          % Bottom
                + ( AP(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*U(:,i,j)...                                                                                                                    % Top
                + AM(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*U(:,i,j+1) + ExpV_y(U(:,i,j), U(:,i,j+1), yc(i,j), yc(i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j)) )*S(4,i,j) ...         % Top 
                + ( AM(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j)) + ImpV_x(U(:,i,j), U(:,i+1,j), xc(i,j), xc(i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j)) )*S(2,i,j)*dU_old(:,i+1,j)...   % Right delta U
                + ( AM(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j)) - ImpV_x(U(:,i,j), U(:,i-1,j), xc(i,j), xc(i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j)) )*S(1,i,j)*dU_old(:,i-1,j) );   % Left delta U 
    
        alpha = a - b*g(:,:,j+1) ;
        g(:,:,j) = alpha\c;
        v(:,1,j) = alpha\(f-b*v(:,1,j+1));
    end
    
    % Bottom 
    

    [UgiB, UgvB, EiB, EvB] = Viscous_Wall_Ghost(U(:,i,1), Face_Normal(1,3,i,1), Face_Normal(2,3,i,1));  


    b = ( AM(U(:,i,1), U(:,i,2), Face_Normal(1,4,i,1), Face_Normal(2,4,i,1)) + ImpV_y(U(:,i,1), U(:,i,2), yc(i,1), yc(i,2), Face_Normal(1,4,i,1), Face_Normal(2,4,i,1)) )*S(4,i,1);     % Top delta U
    
    a = V(i,1)/dt*eye(4) ... 
                + ( AP(U(:,i,1), U(:,i-1,1), Face_Normal(1,1,i,1), Face_Normal(2,1,i,1)) + ImpV_x(U(:,i,1), U(:,i-1,1), xc(i,1), xc(i-1,1), Face_Normal(1,1,i,1), Face_Normal(2,1,i,1)) )*S(1,i,1) ...      % Left
                + ( AP(U(:,i,1), U(:,i+1,1), Face_Normal(1,2,i,1), Face_Normal(2,2,i,1)) - ImpV_x(U(:,i,1), U(:,i+1,1), xc(i,1), xc(i+1,1), Face_Normal(1,2,i,1), Face_Normal(2,2,i,1)) )*S(2,i,1) ...      % Right
                + ( AP(U(:,i,1), UgiB, Face_Normal(1,3,i,1), Face_Normal(2,3,i,1))...                                                                                                                       % Bottom (Wall BC)
                + EiB*AM(U(:,i,1), UgiB, Face_Normal(1,3,i,1), Face_Normal(2,3,i,1))...                                                                                                                    % Bottom (Wall BC)
                + ImpV_by(U(:,i,1), UgvB, EvB, yc(i,1), yc(i,2), Face_Normal(1,3,i,Ny), Face_Normal(2,3,i,Ny)) )*S(3,i,1) ...                                                                             % Bottom (Wall BC)
                + ( AP(U(:,i,1), U(:,i,2), Face_Normal(1,4,i,1), Face_Normal(2,4,i,1)) - ImpV_y(U(:,i,1), U(:,i,2), yc(i,1), yc(i,2), Face_Normal(1,4,i,1), Face_Normal(2,4,i,1)) )*S(4,i,1);               % Top       
    
    f = - (   ( AP(U(:,i,1), U(:,i-1,1), Face_Normal(1,1,i,1), Face_Normal(2,1,i,1))*U(:,i,1)...                                                                                                                    % Left
            + AM(U(:,i,1), U(:,i-1,1), Face_Normal(1,1,i,1), Face_Normal(2,1,i,1))*U(:,i-1,1) + ExpV_x(U(:,i,1), U(:,i-1,1), xc(i,1), xc(i-1,1), Face_Normal(1,1,i,1), Face_Normal(2,1,i,1)) )*S(1,i,1)...          % Left 
            + ( AP(U(:,i,1), U(:,i+1,1), Face_Normal(1,2,i,1), Face_Normal(2,2,i,1))*U(:,i,1)...                                                                                                                    % Right
            + AM(U(:,i,1), U(:,i+1,1), Face_Normal(1,2,i,1), Face_Normal(2,2,i,1))*U(:,i+1,1) + ExpV_x(U(:,i,1), U(:,i+1,1), xc(i,1), xc(i+1,1), Face_Normal(1,2,i,1), Face_Normal(2,2,i,1)) )*S(2,i,1)...          % Right
            + ( AP(U(:,i,1), UgiB, Face_Normal(1,3,i,1), Face_Normal(2,3,i,1))*U(:,i,1)...                                                                                                                          % Bottom (Wall BC)
            + AM(U(:,i,1), UgiB, Face_Normal(1,3,i,1), Face_Normal(2,3,i,1))*UgiB + ExpV_y(U(:,i,1), UgvB, yc(i,1), yc(i,2), Face_Normal(1,3,i,1), Face_Normal(2,3,i,1)) )*S(3,i,1)...                             % Bottom (Wall BC)
            + ( AP(U(:,i,1), U(:,i,2), Face_Normal(1,4,i,1), Face_Normal(2,4,i,1))*U(:,i,1)...                                                                                                                      % Top
            + AM(U(:,i,1), U(:,i,2), Face_Normal(1,4,i,1), Face_Normal(2,4,i,1))*U(:,i,2) + ExpV_y(U(:,i,1), U(:,i,2), yc(i,1), yc(i,2), Face_Normal(1,4,i,1), Face_Normal(2,4,i,1)) )*S(4,i,1)...                  % Top (time n)  
            + ( AM(U(:,i,1), U(:,i+1,1), Face_Normal(1,2,i,1), Face_Normal(2,2,i,1)) + ImpV_x(U(:,i,1), U(:,i+1,1), xc(i,1), xc(i+1,1), Face_Normal(1,2,i,1), Face_Normal(2,2,i,1)) )*S(2,i,1)*dU_old(:,i+1,1)...   % Right delta U
            + ( AM(U(:,i,1), U(:,i-1,1), Face_Normal(1,1,i,1), Face_Normal(2,1,i,1)) - ImpV_x(U(:,i,1), U(:,i-1,1), xc(i,1), xc(i-1,j), Face_Normal(1,1,i,1), Face_Normal(2,1,i,1)) )*S(1,i,1)*dU_old(:,i-1,1) );   % Left delta U 
    
    
    alpha = a - b*g(:,:,2);                
    v(:,1,1) = alpha\(f-b*v(:,1,2));

    dU_new(:,1) = v(:,1,1);
    
    for j = 2:Ny
        dU_new(:,j) = v(:,1,j) - g(:,:,j)*dU_new(:,j-1);
    end

end



%% Function that relaxes the right-most line
function[dU] = Right_Boundary_Tridiagonal_Solver(dU_old, U, U_inlet, Ny, Nx, S, V, Face_Normal, dt, xc, yc)
   
    g = zeros(4,4,Ny); 
    v = zeros(4,1,Ny);
    
    % Top
    [UgiT, UgvT, EiT, EvT] = Outlet(U(:,Nx,Ny), Face_Normal(1,4,Nx,Ny), Face_Normal(2,4,Nx,Ny));

    % [UgiT, UgvT, EiT, EvT] = Inlet(U_inlet, Face_Normal(1,4,Nx,Ny), Face_Normal(2,4,Nx,Ny));
    [UgiR, UgvR, EiR, EvR] = Outlet(U(:,Nx,Ny), Face_Normal(1,2,Nx,Ny), Face_Normal(2,2,Nx,Ny));
 
    a = V(Nx,Ny)/dt*eye(4)...
            + ( AP(U(:,Nx,Ny), U(:,Nx-1,Ny), Face_Normal(1,1,Nx,Ny), Face_Normal(2,1,Nx,Ny)) + ImpV_x(U(:,Nx,Ny), U(:,Nx-1,Ny), xc(Nx,Ny), xc(Nx-1,Ny), Face_Normal(1,1,Nx,Ny), Face_Normal(2,1,Nx,Ny)) )*S(1,Nx,Ny) ...        % Left
            + ( AP(U(:,Nx,Ny), UgiR, Face_Normal(1,2,Nx,Ny), Face_Normal(2,2,Nx,Ny))...                                                                                                                                          % Right (Outflow BC)
            + EiR*AM(U(:,Nx,Ny), UgiR, Face_Normal(1,2,Nx,Ny), Face_Normal(2,2,Nx,Ny))...                                                                                                                                        % Right (Outflow BC)
            + ImpV_bx(U(:,Nx,Ny), UgvR, EvR, xc(Nx,Ny), xc(Nx-1,Ny), Face_Normal(1,2,Nx,Ny), Face_Normal(2,2,Nx,Ny)) )*S(2,Nx,Ny)...                                                                                             % Right (Outflow BC)
            + ( AP(U(:,Nx,Ny), U(:,Nx,Ny-1), Face_Normal(1,3,Nx,Ny), Face_Normal(2,3,Nx,Ny)) + ImpV_y(U(:,Nx,Ny), U(:,Nx,Ny-1), yc(Nx,Ny), yc(Nx,Ny-1), Face_Normal(1,3,Nx,Ny), Face_Normal(2,3,Nx,Ny)) )*S(3,Nx,Ny) ...        % Bottom
            + ( AP(U(:,Nx,Ny), UgiT, Face_Normal(1,4,Nx,Ny), Face_Normal(2,4,Nx,Ny))...                                                                                                                                         % Top (Wall BC)
            + EiT*AM(U(:,Nx,Ny), UgiT, Face_Normal(1,4,Nx,Ny), Face_Normal(2,4,Nx,Ny))...                                                                                                                                      % Top (Wall BC)
            + ImpV_by(U(:,Nx,Ny), UgvT, EvT, yc(Nx,Ny), yc(Nx,Ny-1), Face_Normal(1,4,Nx,Ny), Face_Normal(2,4,Nx,Ny)) )*S(4,Nx,Ny);                                                                                            % Top (Wall BC)
    
    c =  ( AM(U(:,Nx,Ny), U(:,Nx,Ny-1), Face_Normal(1,3,Nx,Ny), Face_Normal(2,3,Nx,Ny)) - ImpV_y(U(:,Nx,Ny), U(:,Nx,Ny-1), yc(Nx,Ny), yc(Nx,Ny-1), Face_Normal(1,3,Nx,Ny), Face_Normal(2,3,Nx,Ny)) )*S(3,Nx,Ny);    % Bottom delta U                                                   
    
    f =  -(   ( AP(U(:,Nx,Ny), U(:,Nx-1,Ny), Face_Normal(1,1,Nx,Ny), Face_Normal(2,1,Nx,Ny))*U(:,Nx,Ny)...                                                                                                                                  % Left
            + AM(U(:,Nx,Ny), U(:,Nx-1,Ny), Face_Normal(1,1,Nx,Ny), Face_Normal(2,1,Nx,Ny))*U(:,Nx-1,Ny) + ExpV_x(U(:,Nx,Ny), U(:,Nx-1,Ny), xc(Nx,Ny), xc(Nx-1,Ny), Face_Normal(1,1,Nx,Ny), Face_Normal(2,1,Nx,Ny)) )*S(1,Nx,Ny)...          % Left
            + ( AP(U(:,Nx,Ny), UgiR, Face_Normal(1,2,Nx,Ny), Face_Normal(2,2,Nx,Ny))*U(:,Nx,Ny)...                                                                                                                                           % Right (Outlet BC)
            + AM(U(:,Nx,Ny), UgiR, Face_Normal(1,2,Nx,Ny), Face_Normal(2,2,Nx,Ny))*UgiR + ExpV_x(U(:,Nx,Ny), UgvR, xc(Nx,Ny), xc(Nx-1,Ny), Face_Normal(1,2,Nx,Ny), Face_Normal(2,2,Nx,Ny)) )*S(2,Nx,Ny)...                                     % Right (Outlet BC)
            + ( AP(U(:,Nx,Ny), U(:,Nx,Ny-1), Face_Normal(1,3,Nx,Ny), Face_Normal(2,3,Nx,Ny))*U(:,Nx,Ny)...                                                                                                                                  % Bottom
            + AM(U(:,Nx,Ny), U(:,Nx,Ny-1), Face_Normal(1,3,Nx,Ny), Face_Normal(2,3,Nx,Ny))*U(:,Nx,Ny-1) + ExpV_y(U(:,Nx,Ny), U(:,Nx,Ny-1), yc(Nx,Ny), yc(Nx,Ny-1), Face_Normal(1,3,Nx,Ny), Face_Normal(2,3,Nx,Ny)) )*S(3,Nx,Ny)...          % Bottom
            + ( AP(U(:,Nx,Ny), UgiT, Face_Normal(1,4,Nx,Ny), Face_Normal(2,4,Nx,Ny))*U(:,Nx,Ny)...                                                                                                                                          % Top (Wall BC)
            + AM(U(:,Nx,Ny), UgiT, Face_Normal(1,4,Nx,Ny), Face_Normal(2,4,Nx,Ny))*UgiT + ExpV_y(U(:,Nx,Ny), UgvT, yc(Nx,Ny), yc(Nx,Ny-1), Face_Normal(1,4,Nx,Ny), Face_Normal(2,4,Nx,Ny)) )*S(4,Nx,Ny) ...                                % Top (Wall BC)  
            + ( AM(U(:,Nx,Ny), U(:,Nx-1,Ny), Face_Normal(1,1,Nx,Ny), Face_Normal(2,1,Nx,Ny)) - ImpV_x(U(:,Nx,Ny), U(:,Nx-1,Ny), xc(Nx,Ny), xc(Nx-1,Ny), Face_Normal(1,1,Nx,Ny), Face_Normal(2,1,Nx,Ny)))*S(1,Nx,Ny)*dU_old(:,Nx-1,Ny) );    % Left delta U  
    

    alpha = a;
    g(:,:,Ny) = alpha\c;
    v(:,1,Ny) = alpha\f;
    
    % Middle

    for j = (Ny-1):-1:2    
        
        [UgiR, UgvR, EiR, EvR] = Outlet(U(:,Nx,j), Face_Normal(1,2,Nx,j), Face_Normal(2,2,Nx,j));

        b = ( AM( U(:,Nx,j), U(:,Nx,j+1), Face_Normal(1,4,Nx,j), Face_Normal(2,4,Nx,j)) + ImpV_y(U(:,Nx,j), U(:,Nx,j+1), yc(Nx,j), yc(Nx,j+1), Face_Normal(1,4,Nx,j), Face_Normal(2,4,Nx,j)) )*S(4,Nx,j);       % Top delta U
    
        a = V(Nx,j)/dt*eye(4) ...
                + ( AP(U(:,Nx,j), U(:,Nx-1,j), Face_Normal(1,1,Nx,j), Face_Normal(2,1,Nx,j)) + ImpV_x(U(:,Nx,j), U(:,Nx-1,j), xc(Nx,j), xc(Nx-1,j), Face_Normal(1,1,Nx,j), Face_Normal(2,1,Nx,j)) )*S(1,Nx,j) ...       % Left
                + ( AP(U(:,Nx,j), UgiR, Face_Normal(1,2,Nx,j), Face_Normal(2,2,Nx,j))...                                                                                                                                 % Right (Outflow BC)
                + EiR*AM(U(:,Nx,j), UgiR, Face_Normal(1,2,Nx,j), Face_Normal(2,2,Nx,j))...                                                                                                                               % Right (Outflow BC)
                + ImpV_bx(U(:,Nx,j), UgvR, EvR, xc(Nx,j), xc(Nx-1,j), Face_Normal(1,2,Nx,j), Face_Normal(2,2,Nx,j)) )*S(2,Nx,j)...                                                                                       % Right (Outflow BC)
                + ( AP(U(:,Nx,j), U(:,Nx,j-1), Face_Normal(1,3,Nx,j), Face_Normal(2,3,Nx,j)) + ImpV_y(U(:,Nx,j), U(:,Nx,j-1), yc(Nx,j), yc(Nx,j-1), Face_Normal(1,3,Nx,j), Face_Normal(2,3,Nx,j)) )*S(3,Nx,j) ...       % Bottom
                + ( AP(U(:,Nx,j), U(:,Nx,j+1), Face_Normal(1,4,Nx,j), Face_Normal(2,4,Nx,j)) - ImpV_y(U(:,Nx,j), U(:,Nx,j+1), yc(Nx,j), yc(Nx,j+1), Face_Normal(1,4,Nx,j), Face_Normal(2,4,Nx,j)) )*S(4,Nx,j);          % Top
    
        c = ( AM(U(:,Nx,j), U(:,Nx,j-1), Face_Normal(1,3,Nx,j), Face_Normal(2,3,Nx,j)) - ImpV_y(U(:,Nx,j), U(:,Nx,j-1), yc(Nx,j), yc(Nx,j-1), Face_Normal(1,3,Nx,j), Face_Normal(2,3,Nx,j)) )*S(3,Nx,j);        % Bottom delta U
    
        f =  -(   ( AP(U(:,Nx,j), U(:,Nx-1,j), Face_Normal(1,1,Nx,j), Face_Normal(2,1,Nx,j))*U(:,Nx,j)...                                                                                                                           $ Left
                + AM(U(:,Nx,j), U(:,Nx-1,j), Face_Normal(1,1,Nx,j), Face_Normal(2,1,Nx,j))*U(:,Nx-1,j) + ExpV_x(U(:,Nx,j), U(:,Nx-1,j), xc(Nx,j), xc(Nx-1,j), Face_Normal(1,1,Nx,j), Face_Normal(2,1,Nx,j)) )*S(1,Nx,j)...          % Left     
                + ( AP(U(:,Nx,j), UgiR, Face_Normal(1,2,Nx,j), Face_Normal(2,2,Nx,j))*U(:,Nx,j)...                                                                                                                                   % Right (Outlet BC)
                + AM(U(:,Nx,j), UgiR, Face_Normal(1,2,Nx,j), Face_Normal(2,2,Nx,j))*UgiR + ExpV_x(U(:,Nx,j), UgvR, xc(Nx,j), xc(Nx-1,j), Face_Normal(1,2,Nx,j), Face_Normal(2,2,Nx,j)) )*S(2,Nx,j)...                                  % Right (Outlet BC)                
                + ( AP(U(:,Nx,j), U(:,Nx,j-1), Face_Normal(1,3,Nx,j), Face_Normal(2,3,Nx,j))*U(:,Nx,j)...                                                                                                                           % Bottom
                + AM(U(:,Nx,j), U(:,Nx,j-1), Face_Normal(1,3,Nx,j), Face_Normal(2,3,Nx,j))*U(:,Nx,j-1) + ExpV_y(U(:,Nx,j), U(:,Nx,j-1), yc(Nx,j), yc(Nx,j-1), Face_Normal(1,3,Nx,j), Face_Normal(2,3,Nx,j)) )*S(3,Nx,j)...          % Bottom
                + ( AP(U(:,Nx,j), U(:,Nx,j+1), Face_Normal(1,4,Nx,j), Face_Normal(2,4,Nx,j))*U(:,Nx,j)...                                                                                                                           % Top
                + AM(U(:,Nx,j), U(:,Nx,j+1), Face_Normal(1,4,Nx,j), Face_Normal(2,4,Nx,j))*U(:,Nx,j+1) + ExpV_y(U(:,Nx,j), U(:,Nx,j+1), yc(Nx,j), yc(Nx,j+1), Face_Normal(1,4,Nx,j), Face_Normal(2,4,Nx,j)) )*S(4,Nx,j) ...         % Top
                + ( AM(U(:,Nx,j), U(:,Nx-1,j), Face_Normal(1,1,Nx,j), Face_Normal(2,1,Nx,j)) - ImpV_x(U(:,Nx,j), U(:,Nx-1,j), xc(Nx,j), xc(Nx-1,j), Face_Normal(1,1,Nx,j), Face_Normal(2,1,Nx,j)) )*S(1,Nx,j)*dU_old(:,Nx-1,j) );   % Left delta U
    
        alpha = a - b*g(:,:,j+1);
        g(:,:,j) = alpha\c;
        v(:,1,j) = alpha\(f-b*v(:,1,j+1));
    
    end
    
    % Bottom

    [UgiB, UgvB, EiB, EvB] = Viscous_Wall_Ghost(U(:,Nx,1), Face_Normal(1,3,Nx,1), Face_Normal(2,3,Nx,1));
    [UgiR, UgvR, EiR, EvR] = Outlet(U(:,Nx,1), Face_Normal(1,2,Nx,1), Face_Normal(2,2,Nx,1));
    
    b = ( AM(U(:,Nx,1), U(:,Nx,2), Face_Normal(1,4,Nx,1), Face_Normal(2,4,Nx,1)) + ImpV_y( U(:,Nx,1), U(:,Nx,2), yc(Nx,1), yc(Nx,2), Face_Normal(1,4,Nx,1), Face_Normal(2,4,Nx,1)) )*S(4,Nx,1);           % Top delta U
    
    a = V(Nx,1)/dt*eye(4) ...
            + ( AP(U(:,Nx,1), U(:,Nx-1,1), Face_Normal(1,1,Nx,1), Face_Normal(2,1,Nx,1)) + ImpV_x(U(:,Nx,1), U(:,Nx-1,1), xc(Nx,1), xc(Nx-1,1), Face_Normal(1,1,Nx,1), Face_Normal(2,1,Nx,1)) )*S(1,Nx,1) ...       % Left
            + ( AP(U(:,Nx,1), UgiR, Face_Normal(1,2,Nx,1), Face_Normal(2,2,Nx,1))...                                                                                                                                 % Right (Outflow BC)
            + EiR*AM(U(:,Nx,1), UgiR, Face_Normal(1,2,Nx,1), Face_Normal(2,2,Nx,1))...                                                                                                                               % Right (Ourflow BC)
            + ImpV_bx(U(:,Nx,1), UgvR, EvR, xc(Nx,1), xc(Nx-1,1), Face_Normal(1,2,Nx,1), Face_Normal(2,2,Nx,1)) )*S(2,Nx,1)...                                                                                       % Right (Outflow BC)
            + ( AP(U(:,Nx,1), UgiB, Face_Normal(1,3,Nx,1), Face_Normal(2,3,Nx,1))...                                                                                                                                % Bottom (Wall BC)
            + EiB*AM(U(:,Nx,1), UgiB, Face_Normal(1,3,Nx,1), Face_Normal(2,3,Nx,1))...                                                                                                                             % Bottom (Wall BC)
            + ImpV_by(U(:,Nx,1), UgvB, EvB, yc(Nx,1), yc(Nx,2), Face_Normal(1,3,Nx,1), Face_Normal(2,3,Nx,1)) )*S(3,Nx,1) ...                                                                                     % Bottom (Wall BC) 
            + ( AP(U(:,Nx,1), U(:,Nx,2), Face_Normal(1,4,Nx,1), Face_Normal(2,4,Nx,1)) - ImpV_y(U(:,Nx,1), U(:,Nx,2), yc(Nx,1), yc(Nx,2), Face_Normal(1,4,Nx,1), Face_Normal(2,4,Nx,1)) )*S(4,Nx,1);                % Top 

    f =  -(   ( AP(U(:,Nx,1), U(:,Nx-1,1), Face_Normal(1,1,Nx,1), Face_Normal(2,1,Nx,1))*U(:,Nx,1)...                                                                                                                           % Left
            + AM(U(:,Nx,1), U(:,Nx-1,1), Face_Normal(1,1,Nx,1), Face_Normal(2,1,Nx,1))*U(:,Nx-1,1) + ExpV_x(U(:,Nx,1), U(:,Nx-1,1), xc(Nx,1), xc(Nx-1,1), Face_Normal(1,1,Nx,1), Face_Normal(2,1,Nx,1)) )*S(1,Nx,1)...          % Left                        
            + ( AP(U(:,Nx,1), UgiR, Face_Normal(1,2,Nx,1), Face_Normal(2,2,Nx,1))*U(:,Nx,1)...                                                                                                                                   % Right (Outlet BC)
            + AM(U(:,Nx,1), UgiR, Face_Normal(1,2,Nx,1), Face_Normal(2,2,Nx,1))*UgiR + ExpV_x(U(:,Nx,1), UgvR, xc(Nx,1), xc(Nx-1,1), Face_Normal(1,2,Nx,1), Face_Normal(2,2,Nx,1)) )*S(2,Nx,1)...                            % Right (Outlet BC)
            + ( AP(U(:,Nx,1), UgiB, Face_Normal(1,3,Nx,1), Face_Normal(2,3,Nx,1))*U(:,Nx,1)...                                                                                                                                  % Bottom (Wall BC)
            + AM(U(:,Nx,1), UgiB, Face_Normal(1,3,Nx,1), Face_Normal(2,3,Nx,1))*UgiB + ExpV_y(U(:,Nx,1), UgvB, yc(Nx,1), yc(Nx,2), Face_Normal(1,3,Nx,1), Face_Normal(2,3,Nx,1)) )*S(3,Nx,1)...                                % Bottom (Wall BC)
            + ( AP(U(:,Nx,1), U(:,Nx,2), Face_Normal(1,4,Nx,1), Face_Normal(2,4,Nx,1))*U(:,Nx,1)...                                                                                                                             % Top
            + AM(U(:,Nx,1), U(:,Nx,2), Face_Normal(1,4,Nx,1), Face_Normal(2,4,Nx,1))*U(:,Nx,2) + ExpV_y(U(:,Nx,1), U(:,Nx,2), yc(Nx,1), yc(Nx,2), Face_Normal(1,4,Nx,1), Face_Normal(2,4,Nx,1)) )*S(4,Nx,1)...                  % Top              
            + ( AM(U(:,Nx,1), U(:,Nx-1,1), Face_Normal(1,1,Nx,1), Face_Normal(2,1,Nx,1)) - ImpV_x(U(:,Nx,1), U(:,Nx-1,1), xc(Nx,1), xc(Nx-1,1), Face_Normal(1,1,Nx,1), Face_Normal(2,1,Nx,1)) )*S(1,Nx,1)*dU_old(:,Nx-1,1) );   % Left delta U
    
    
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
%% Function that calculates U for ghost cells at a viscous wall
function [Ugiw, Ugvw, Ei, Ev] = Viscous_Wall_Ghost(U, nx, ny)
    
    Inviscid_Ghost_primitive = zeros(4,1);
    Inside_Primitives = consToPrim(U);     
    Q = Viscous_Primitives(U);
    
    rho = Inside_Primitives(1);
    u = Inside_Primitives(2);
    v = Inside_Primitives(3);
    p = Inside_Primitives(4);

    Ti = Q(4);

    Inviscid_Ghost_primitive(1) = rho;
    Inviscid_Ghost_primitive(2) = u - 2*(u*nx + v*ny)*nx;
    Inviscid_Ghost_primitive(3) = v - 2*(u*nx + v*ny)*ny;
    Inviscid_Ghost_primitive(4) = p;
    
    Tw = 300;

    rho_visc = 2*(rho*Ti/Tw)-rho;
    u_visc = -u;
    v_visc = -v;
    p_visc = p;

    Viscid_Ghost_primitive(1) = rho_visc;
    Viscid_Ghost_primitive(2) = u_visc;
    Viscid_Ghost_primitive(3) = v_visc;
    Viscid_Ghost_primitive(4) = p_visc;
    
    Ugiw = primToCons(Inviscid_Ghost_primitive);
    Ugvw = primToCons(Viscid_Ghost_primitive);

    Ei = [1, 0, 0, 0;
         0, (1 - 2*nx^2), -2*nx*ny, 0;
         0, -2*nx*ny,  (1-2*ny^2), 0;
         0, 0, 0, 1];    

    Ev = [2*Ti/Tw-1, 0, 0, 2*rho/Tw;
        0, -1, 0, 0;
        0, 0, -1, 0;
        0, 0, 0, -1];

end

%% Function that makes symmetry ghost cell
function [Ugiw, Ugvw, Ei, Ev] = symmetry(U, nx, ny)
    
    Inviscid_Ghost_primitive = zeros(4,1);
    Inside_Primitives = consToPrim(U);     
    Q = Viscous_Primitives(U);
    
    rho = Inside_Primitives(1);
    u = Inside_Primitives(2);
    v = Inside_Primitives(3);
    p = Inside_Primitives(4);

    Inviscid_Ghost_primitive(1) = rho;
    Inviscid_Ghost_primitive(2) = u - 2*(u*nx + v*ny)*nx;
    Inviscid_Ghost_primitive(3) = v - 2*(u*nx + v*ny)*ny;
    Inviscid_Ghost_primitive(4) = p;

    
    Ugiw = primToCons(Inviscid_Ghost_primitive);
    Ugvw = primToCons(Inviscid_Ghost_primitive);

    Ei = [1, 0, 0, 0;
         0, (1 - 2*nx^2), -2*nx*ny, 0;
         0, -2*nx*ny,  (1-2*ny^2), 0;
         0, 0, 0, 1];    

    Ev = [1, 0, 0, 0;
         0, (1 - 2*nx^2), -2*nx*ny, 0;
         0, -2*nx*ny,  (1-2*ny^2), 0;
         0, 0, 0, 1];    


end

%% Function that calculates outlet boundary conditions
function [Ugiw, Ugvw, Ei,Ev] = Outlet(U, nx, ny)
        
    Ugiw = U;
    Ugvw = U;

    Ei = eye(4); 
    Ev = eye(4);

end
%% Function that calculate inlet boundary conditions
function [Ugii, Ugvi, Eii,Evi] = Inlet(U_inlet, nx, ny)
        
    Ugii = U_inlet;
    Ugvi = U_inlet;

    Eii = zeros(4,4); 
    Evi = zeros(4,4);

end

%% Function that calculates residual for GSLR
function [residual] = calculateLocalResidual(dU, dU_old, U, V, dt, Ny, S, Face_Normal, i, xc, yc)
    
    lineresidual = 0;


    for j = 2:Ny-1

        b = ( AM(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j)) + ImpV_y(U(:,i,j+1), U(:,i,j), yc(i,j), yc(i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j)) )*S(4,i,j); 

        a = V(i,j)/dt*eye(4) ...
                    + ( AP(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j)) + ImpV_x(U(:,i,j), U(:,i-1,j), xc(i,j), xc(i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j)) )*S(1,i,j) ... % Left
                    + ( AP(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j)) - ImpV_x(U(:,i,j), U(:,i+1,j), xc(i,j), xc(i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j)) )*S(2,i,j) ... % Right
                    + ( AP(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j)) + ImpV_y(U(:,i,j), U(:,i,j-1), yc(i,j), yc(i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j)) )*S(3,i,j) ... % Bottom
                    + ( AP(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j)) - ImpV_y(U(:,i,j), U(:,i,j+1), yc(i,j), yc(i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j)) )*S(4,i,j);    % Top 
    
        c = ( AM(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j)) - ImpV_y(U(:,i,j-1), U(:,i,j), yc(i,j), yc(i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j)) )*S(3,i,j); 
    
        f =  -( ( AP(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*U(:,i,j)...
                + AM(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*U(:,i-1,j) + ExpV_x(U(:,i,j), U(:,i-1,j), xc(i,j), xc(i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j)) )*S(1,i,j)...   % Left 
                + ( AP(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*U(:,i,j)...
                + AM(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*U(:,i+1,j) + ExpV_x(U(:,i,j), U(:,i+1,j), xc(i,j), xc(i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j)) )*S(2,i,j)...   % Right 
                + ( AP(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*U(:,i,j)...
                + AM(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*U(:,i,j-1) + ExpV_y(U(:,i,j), U(:,i,j-1), yc(i,j), yc(i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j)) )*S(3,i,j)...   % Bottom 
                + ( AP(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*U(:,i,j)...
                + AM(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*U(:,i,j+1) + ExpV_y(U(:,i,j), U(:,i,j+1), yc(i,j), yc(i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j)) )*S(4,i,j) ...  % Top              
                + ( AM(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j)) + ImpV_x(U(:,i+1,j), U(:,i,j), xc(i,j), xc(i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j)) )*S(2,i,j)*dU_old(:,i+1,j)...
                + ( AM(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j)) - ImpV_x(U(:,i-1,j), U(:,i,j), xc(i,j), xc(i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j)) )*S(1,i,j)*dU_old(:,i-1,j) );      % Left (time n) and right (time n-1) delta U 
    
                res =  b*dU(:,i,j+1) + a*dU(:,i,j) + c*dU(:,i,j-1) - f;      
                lineresidual = lineresidual + res(1)^2;

    end

    residual = sqrt(lineresidual);

end

%% Function that calculates outer residual
function [total_residual] = Calculate_Residual(U, V, S, Face_Normal, Nx, Ny, xc, yc)

    residual = 0;
    
    for i=2:Nx-1
        for j = 2:Ny-1
    
            residual =  residual + ((( AP(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*U(:,i,j)...
                + AM(U(:,i,j), U(:,i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j))*U(:,i-1,j) + ExpV_x(U(:,i,j), U(:,i-1,j), xc(i,j), xc(i-1,j), Face_Normal(1,1,i,j), Face_Normal(2,1,i,j)) )*S(1,i,j)...   
                + ( AP(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*U(:,i,j)...
                + AM(U(:,i,j), U(:,i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j))*U(:,i+1,j) + ExpV_x(U(:,i,j), U(:,i+1,j), xc(i,j), xc(i+1,j), Face_Normal(1,2,i,j), Face_Normal(2,2,i,j)) )*S(2,i,j)...   
                + ( AP(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*U(:,i,j)...
                + AM(U(:,i,j), U(:,i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j))*U(:,i,j-1) + ExpV_y(U(:,i,j), U(:,i,j-1), yc(i,j), yc(i,j-1), Face_Normal(1,3,i,j), Face_Normal(2,3,i,j)) )*S(3,i,j)...  
                + ( AP(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*U(:,i,j)...
                + AM(U(:,i,j), U(:,i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j))*U(:,i,j+1) + ExpV_y(U(:,i,j), U(:,i,j+1), yc(i,j), yc(i,j+1), Face_Normal(1,4,i,j), Face_Normal(2,4,i,j)) )*S(4,i,j) )/V(i,j)).^2;      % Top (time n)               
               
        end
    end
   
    res = residual(1);

    total_residual = sqrt(res);

end

%% Function that calculates viscous primitives (rho, u, v, T)
function [L] = Viscous_Primitives(U)
    
    R = 287;
    gamma = 1.4;
    rho = U(1);          % Density
    rhou = U(2);         % rho * u
    rhov = U(3);         % rho * v
    E = U(4);            % Total energy

    % Calculate primitive variables
    L = zeros(4,1);
    L(1) = rho;                  % Density is directly available
    L(2) = rhou / rho;          % Velocity in x-direction: u = (rho * u) / rho
    L(3) = rhov / rho;          % Velocity in y-direction: v = (rho * v) / rho

    % Total energy equation to calculate pressure
    kinetic_energy = 0.5 * rho * (L(2)^2 + L(3)^2);
    P = (gamma - 1) * (E - kinetic_energy);  % Pressure: p = (gamma - 1) * (E - 0.5*rho*(u^2 + v^2))
    L(4) = P/(R*rho);


end

%% Function that calculates implicit viscous fluxes in x direction
function [Fv] = ImpV_x(U, Uii, d1, d2, nx, ny)

    Vi = Viscous_Primitives(U);
    Vii = Viscous_Primitives(Uii);

    cv = 717;
    rho = 0.5*(Vi(1) + Vii(1));
    u = 0.5*(Vi(2) + Vii(2));
    v = 0.5*(Vi(3) + Vii(3));
    T = 0.5*(Vi(4) + Vii(4));
    mu = 1.458*10^(-6)*T^(3/2)/(T+110.3);
    lambda = -2/3*mu;
    cp = 1005;
    Pr = 0.71;
    k = cp*mu/Pr;

    Mxx = [0, 0, 0, 0;
           0, lambda + 2*mu, 0, 0;
           0, 0, mu, 0;
           0, u*(lambda + 2*mu), v*mu, k];
    
    N = [1, 0, 0, 0;
        -u/rho, 1/rho, 0, 0;
        -v/rho, 0, 1/rho, 0;
        (u^2+v^2)/(2*cv*rho) - T/rho, -u/(cv*rho), -v/(cv*rho), 1/(cv*rho)];
    
    
    dn = abs(d2-d1);

    Fv = -(Mxx)*N/dn*nx;

end

%% Function that calculates implicit viscous fluxes in y direction
function [Fv] = ImpV_y(U, Uii, d1, d2, nx, ny)

    Vi = Viscous_Primitives(U);
    Vii = Viscous_Primitives(Uii);

    cv = 717;
    rho = 0.5*(Vi(1) + Vii(1));
    u = 0.5*(Vi(2) + Vii(2));
    v = 0.5*(Vi(3) + Vii(3));
    T = 0.5*(Vi(4) + Vii(4));
    mu = 1.458*10^(-6)*T^(3/2)/(T+110.3);
    lambda = -2/3*mu;
    cp = 1005;
    Pr = 0.71;
    k = cp*mu/Pr;

    Myy = [0, 0, 0, 0;
           0, mu, 0, 0;
           0, 0, lambda + 2*mu, 0;
           0, u*mu, v*(lambda + 2*mu), k];
    
    N = [1, 0, 0, 0;
        -u/rho, 1/rho, 0, 0;
        -v/rho, 0, 1/rho, 0;
        (u^2+v^2)/(2*cv*rho) - T/rho, -u/(cv*rho), -v/(cv*rho), 1/(cv*rho)];
    
    
    dn = abs(d2-d1);

    Fv = -Myy*N/dn*ny;

end

%% Function that calculates implicit viscous boundary fluxes in x direction
function [Fv] = ImpV_bx(U, Uii, E, d1, d2, nx, ny)

    Vi = Viscous_Primitives(U);
    Vii = Viscous_Primitives(Uii);

    cv = 717;
    rho = 0.5*(Vi(1) + Vii(1));
    u = 0.5*(Vi(2) + Vii(2));
    v = 0.5*(Vi(3) + Vii(3));
    T = 0.5*(Vi(4) + Vii(4));
    mu = 1.458*10^(-6)*T^(3/2)/(T+110.3);
    lambda = -2/3*mu;
    cp = 1005;
    Pr = 0.71;
    k = cp*mu/Pr;

    Mxx = [0, 0, 0, 0;
           0, lambda + 2*mu, 0, 0;
           0, 0, mu, 0;
           0, u*(lambda + 2*mu), v*mu, k];
    
    N = [1, 0, 0, 0;
        -u/rho, 1/rho, 0, 0;
        -v/rho, 0, 1/rho, 0;
        (u^2+v^2)/(2*cv*rho) - T/rho, -u/(cv*rho), -v/(cv*rho), 1/(cv*rho)];
    
    
    dn = abs(d2-d1);

    Fv = -Mxx/dn*(E-eye(4))*N*nx;

end

%% Function that calculates implicit boundary fluxes in y direction
function [Fv] = ImpV_by(U, Uii, E, d1, d2, nx, ny)

    Vi = Viscous_Primitives(U);
    Vii = Viscous_Primitives(Uii);

    cv = 717;
    rho = 0.5*(Vi(1) + Vii(1));
    u = 0.5*(Vi(2) + Vii(2));
    v = 0.5*(Vi(3) + Vii(3));
    T = 0.5*(Vi(4) + Vii(4));
    mu = 1.458*10^(-6)*T^(3/2)/(T+110.3);
    lambda = -2/3*mu;
    cp = 1005;
    Pr = 0.71;
    k = cp*mu/Pr;

    Myy = [0, 0, 0, 0;
           0, mu, 0, 0;
           0, 0, lambda + 2*mu, 0;
           0, u*mu, v*(lambda + 2*mu), k];
    
    N = [1, 0, 0, 0;
        -u/rho, 1/rho, 0, 0;
        -v/rho, 0, 1/rho, 0;
        (u^2+v^2)/(2*cv*rho) - T/rho, -u/(cv*rho), -v/(cv*rho), 1/(cv*rho)];
    
    
    dn = abs(d2-d1);

    Fv = -Myy/dn*(E-eye(4))*N*ny;

end

%% Function that calculates explicit viscous fluxes in y direction
function [Fv] = ExpV_y(Ui, Uii, d1, d2, nx, ny)
    
    Vi = Viscous_Primitives(Ui);
    Vii = Viscous_Primitives(Uii);
    
    u = 0.5*(Vi(2) + Vii(2));
    v = 0.5*(Vi(3) + Vii(3));
    T = 0.5*(Vi(4) + Vii(4));

    mu = 1.458*10^(-6)*T^(3/2)/(T+110.3);
    lambda = -2/3*mu;
    cp = 1005;
    Pr = 0.71;
    k = cp*mu/Pr;
    

    Myy = [0, 0, 0, 0;
           0, mu, 0, 0;
           0, 0, lambda + 2*mu, 0;
           0, u*mu, v*(lambda + 2*mu), k];
    
    
    dn = abs(d2-d1);
    dvdn = (Vii - Vi)/dn*ny;

    Fv = -Myy*dvdn*ny;


end

%% Function that calculates explicit viscous fluxes in x direction
function [Fv] = ExpV_x(Ui, Uii, d1, d2, nx, ny)
    
    Vi = Viscous_Primitives(Ui);
    Vii = Viscous_Primitives(Uii);
    

    u = 0.5*(Vi(2) + Vii(2));
    v = 0.5*(Vi(3) + Vii(3));
    T = 0.5*(Vi(4) + Vii(4));

    mu = 1.458*10^(-6)*T^(3/2)/(T+110.3);
    lambda = -2/3*mu;
    cp = 1005;
    Pr = 0.71;
    k = cp*mu/Pr;

    Mxx = [0, 0, 0, 0;
           0, lambda + 2*mu, 0, 0;
           0, 0, mu, 0;
           0, u*(lambda + 2*mu), v*mu, k];

    
    dn = abs(d2-d1);
    dvdn = (Vii-Vi)/dn*nx;

    Fv = -Mxx*dvdn*nx;

end



%% Calculate dt
function [dt] = Calculate_dt(U, S, Nx, Ny)

    gamma = 1.4;
    int_dt = zeros(Nx,Ny);
    
    for i = 1:Nx
        for j = 1:Ny

            V = consToPrim(U(:,i,j));
            dx = min(min(S(1,i,j), S(2,i,j)));
            dy = min(min(S(3,i,j), S(4,i,j)));
            c = sqrt(gamma*V(4)/V(1));

            int_dt(i,j) =  1/( abs(V(2)/dx) + abs(V(3)/dy) + c*sqrt(1/dx^2 + 1/dy^2) );

        end
    end

    dt = min(min(int_dt));
    
end
