clear all;
close all;
clc;

%% Grid parameters
Ny = 100;
Nx = 200;
Nr = Ny + 1;                % Number of points in radial direction
Ntheta = Nx + 1;            % Number of points revolving around theta
Cylinder_Radius = 0.1;  % Radius of cylinder mesh is fitted around
R1 = 0.3;              % Radius at theta = 0 
R2 = 0.45;              % Radius at theta = pi/2 and 3pi/2
dr_min = 0.001;         % Thickness of first layer
theta_i = 3*pi/2;         % Starting angle for theta sweep  
theta_f = pi/2;       % Ending angle for theta sweep
    
[x,y] = GridGenerator(Nr, Ntheta, Cylinder_Radius, R1, R2, dr_min, theta_i, theta_f);
[center, S, V, Face_Normal] = GridGeometry(x,y,Nr,Ntheta);

clearvars -except center V S Face_Normal x y Nx Ny;


%% Functions

% This function uses the newton method
function [k] = NewtonMethod(max_dist, n_points, d_min)

k = 1;
knew = 1/2;
ratio = abs(knew-k);

    while ratio >= 0.00000000001
        func = d_min - max_dist*(exp(k/(n_points-1)) - 1)/(exp(k) - 1);
        func_prime = -max_dist*( ((1/(n_points-1)*exp(k/(n_points-1)))*(exp(k)-1) - (exp(k/(n_points-1)) - 1)*exp(k) )/(exp(k) - 1)^2 );
        knew = k - func/func_prime;
        ratio = abs(k-knew);
        k = knew;
    end

end

% This function generates the grid and stretches it in the radial direction
% while also smoothly extending the radius from R1 to R2 as we speed from
% theta = 0 to theta = pi/2 or 3pi/2
function [x,y] = GridGenerator(Nr, Ntheta, Cylinder_Radius, R1, R2, dr_min, theta_i, theta_f)
    
    x = zeros(Ntheta,Nr);
    y = zeros(Ntheta,Nr);     
    theta = zeros(Ntheta,1);
    r = zeros(Ntheta, Nr);

    dtheta = (theta_f - theta_i)/(Ntheta-1);
    r(:,1) = Cylinder_Radius;
    
    % Create points sweeping through theta
    for i = 1:Ntheta
        theta(i) = theta_i + (i-1)*dtheta;
    end    

    % Create r points
    for i = 1:Ntheta

        R_max = R1 + (R2 - R1)*cos(theta(i)); % Calculates max R based on angle
        k1 = NewtonMethod(R_max,Nr,dr_min);

        for j=2:Nr
            r(i,j) = r(i,1) + R_max*( (exp(k1*(j-1)/(Nr-1)) - 1)/(exp(k1)-1));    
        end

    end    
    
    % Transfer to x and y coordinates
    for j = 1:Nr
        for i = 1:Ntheta
            x(i,j) = r(i,j)*cos(theta(i));
            y(i,j) = r(i,j)*sin(theta(i));
        end
    end

end

function [center, S, V, Face_Normal] = GridGeometry(x,y, Nr, Ntheta)

    center = zeros(2,Ntheta-1,Nr-1);
    S = zeros(4,Ntheta-1,Nr-1);
    V = zeros(Ntheta-1,Nr-1);
    Face_Normal = zeros(1,3,Ntheta-1,Nr-1);
    
    for i = 1:Ntheta-1
        for j = 1:Nr-1
        
        % x and y center of cell where indices are: center[(1=x, 2=y), column, row]
        center(1,i,j) = (x(i,j) + x(i+1,j) + x(i,j+1) + x(i+1,j+1))/4;
        center(2,i,j) = (y(i,j) + y(i+1,j) + y(i,j+1) + y(i+1,j+1))/4;
        
        % Vectors from one points to the next
        DA = [x(i,j) - x(i+1,j), y(i,j) - y(i+1,j)]; % Bottom side (y varies)
        AB = [x(i,j+1) - x(i,j), y(i,j+1) - y(i,j)]; % Left side (x varies)
        BC = [x(i+1,j+1) - x(i,j+1), y(i+1,j+1) - y(i,j+1)]; %  Top side (y varies)
        CD = [x(i+1,j) - x(i+1,j+1), y(i+1,j) - y(i+1,j+1)]; %  Right side (x varies)    
        
        % Face area calculatation where indices are: S[(1=Left, 2=Bottom, 3=Right, 4=Top), column, row]
        S(3,i,j) = sqrt(DA(1)^2 + DA(2)^2); % Bottom side area
        S(4,i,j) = sqrt(BC(1)^2 + BC(2)^2); % Top side area
        S(2,i,j) = sqrt(CD(1)^2 + CD(2)^2); % Right side area
        S(1,i,j) = sqrt(AB(1)^2 + AB(2)^2); % Left side area
        
        % Volume calculation using cross product of vectors
        V(i,j) = 1/2*abs(AB(1)*DA(2) - AB(2)*DA(1)) + 1/2*abs(CD(1)*BC(2) - CD(2)*BC(1));

        % Face normal calculation where indices are: Face_Normal[(1=x, 2=y), (1=L, 2=B, 3=R, 4=T), column, row)]
        Face_Normal(1,3,i,j) = -DA(2)/abs(S(3,i,j));
        Face_Normal(1,1,i,j) = -AB(2)/abs(S(1,i,j));
        Face_Normal(1,4,i,j) = -BC(2)/abs(S(4,i,j));
        Face_Normal(1,2,i,j) = -CD(2)/abs(S(2,i,j));
    
        Face_Normal(2,3,i,j) = DA(1)/abs(S(3,i,j));
        Face_Normal(2,1,i,j) = AB(1)/abs(S(1,i,j));
        Face_Normal(2,4,i,j) = BC(1)/abs(S(4,i,j));
        Face_Normal(2,2,i,j) = CD(1)/abs(S(2,i,j));
        
        end
    end

end
