clear all;
close all;
clc;

%% Parameters (Change here!)
L1 = 10;             % Length of section before ramp
L2 = 10;             % Length of ramp section
L3 = 10;             % Length of section after ramp.
Nx = 800;           % Number of cells in x-direction
Ny = 400;            % Number of cells in y-direction
Ramp_Angle = 15;    % Angle of ramp in degrees
inlet_height = 10;   % Height of inlet

% Makes sure the angle and length of the ramp do go higher than the inlet
% height

check = L2*tand(Ramp_Angle);

if check >= inlet_height
    msg = sprintf('Inlet must be larger than %d ', check);   
    disp(msg);
    inlet_height = input('New height: ');    
end

%% Code 
[x,y] = GridGenerator(L1,L2,L3,Nx,Ny,Ramp_Angle,inlet_height);
[center, S, V, Face_Normal] = GridGeometry(x,y, Nx, Ny);

centerx = squeeze(center(1,:,:));
centery = squeeze(center(2,:,:));

clearvars -except center V S Face_Normal x y Nx Ny;

%% Functions

% This function generates the grid as uniformly as possiblem, though at
% boundaries the dx and dy may be adjusted to snap to the stated lenghts.

function [x,y] =  GridGenerator(L1, L2, L3, Nx, Ny, Ramp_Angle, inlet_height)
    
    L_tot = L1 + L2 + L3; % Total Length of grid
    dx = L_tot /(Nx+1);   % Uniform dx
    
    % Initialize x and y matrices
    x = zeros(Nx+1, Ny+1);
    y = zeros(Nx+1, Ny+1);
    
    % Fill in x's
    for i = 1:Nx+1
        x(i,:) = (i-1)*dx;
    end
    
    % Snap values to specific locations so that segments are lengths
    % stated
    for i = 2:Nx+1
        if x(i,1) > L1 && x(i,1) < L1 + dx
            x(i,:) = L1;
        elseif x(i,1) > L1 + L2 && x(i,1) < L1 + L2 + dx
            x(i,:) = L1 + L2;
        elseif i == Nx + 1
            x(i,:) = L_tot;
        end
    end

    %% Generate the grid
    for i = 1:Nx+1
        
        L = x(i,1); % Current length

        % Flat region (first segment)
        if L <= L1       

            dy = inlet_height/(Ny);
            for j = 1:Ny+1
                y(i,j) = (j-1)*dy;
            end   

        % Ramp region (second segment)
        elseif L > L1 && L <= (L1 + L2)             
           
            % This if-statement checks if the current x-value was snapped
            % to L1 or (L1 + L2) and if so it changes the calculation for
            % dy_ramp since it is a function of dx which changes at these
            % boundaries.

            if x(i-1,1) == L1                
                dy_ramp = (x(i,1) - x(i-1,1))*tand(Ramp_Angle);
            elseif L + dx > L1 + L2 && L < L1 + L2 + dx  
                dy_ramp = (x(i,1) - x(i-1,1))*tand(Ramp_Angle);
            else 
                dy_ramp = dx*tand(Ramp_Angle);
            end
            
            y(i,1) = y(i-1,1) + dy_ramp;  % y at j = 1
            dy = (inlet_height - y(i,1))/(Ny);  % Recalculate dy due to ramp height change

            for j = 2:Ny+1
                y(i,j) = y(i,1) + (j-1)*dy;  % Set y for each row after j = 1            
            end
           
        % Flat region (third segment)
        else            
            y(i,1) = L2*tand(Ramp_Angle); %y at j = 1
            dy = (inlet_height - L2*tand(Ramp_Angle))/(Ny); 

            for j = 2:Ny+1
                y(i,j) = L2*tand(Ramp_Angle) + (j-1)*dy;
            end

        end
    end

end

% This function finds the center, the face areas, the volume,
% and the face normals of each cell.
function [center, S, V, Face_Normal] = GridGeometry(x,y, Nx, Ny)

center = zeros(2,Nx,Ny);
S = zeros(4,Nx,Ny);
V = zeros(Nx,Ny);
Face_Normal = zeros(2,4,Nx,Ny);

for i = 1:Nx
    for j = 1:Ny
    
    % x and y center of cell where indices are: center[(1=x, 2=y), column, row]
    center(1,i,j) = (x(i,j) + x(i+1,j) + x(i,j+1) + x(i+1,j+1))/4;
    center(2,i,j) = (y(i,j) + y(i+1,j) + y(i,j+1) + y(i+1,j+1))/4;
    
    % Vectors from one points to the next
    DA = [x(i,j) - x(i,j+1), y(i,j) - y(i,j+1)]; 
    AB = [x(i+1,j) - x(i,j), y(i+1,j) - y(i,j)];
    BC = [x(i+1,j+1) - x(i+1,j), y(i+1,j+1) - y(i+1,j)];
    CD = [x(i,j+1) - x(i+1,j+1), y(i,j+1) - y(i+1,j+1)];
    
    % Face area calculatation where indices are: S[(1=Left, 2=Right, 3=Bottom, 4=Top), column, row]
    S(1,i,j) = sqrt(DA(1)^2 + DA(2)^2); % Left side area
    S(3,i,j) = sqrt(AB(1)^2 + AB(2)^2); % Bottom side area
    S(2,i,j) = sqrt(BC(1)^2 + BC(2)^2); % Right side area
    S(4,i,j) = sqrt(CD(1)^2 + CD(2)^2); % Top side area
    
    % Volume calculation using cross product of vectors
    V(i,j) = 1/2*abs(DA(1)*AB(2) - DA(2)*AB(1)) + 1/2*abs(BC(1)*CD(2) - BC(2)*CD(1));
    
    % Face normal calculation where indices are: Face_Normal[(1=x, 2=y), (1=L, 2=R, 3=B, 4=T), column, row)]
    Face_Normal(1,1,i,j) = DA(2)/abs(S(1,i,j));
    Face_Normal(1,3,i,j) = AB(2)/abs(S(3,i,j));
    Face_Normal(1,2,i,j) = BC(2)/abs(S(2,i,j));
    Face_Normal(1,4,i,j) = CD(2)/abs(S(4,i,j));

    Face_Normal(2,1,i,j) = DA(1)/abs(S(1,i,j));
    Face_Normal(2,3,i,j) = -AB(1)/abs(S(3,i,j));
    Face_Normal(2,2,i,j) = BC(1)/abs(S(2,i,j));
    Face_Normal(2,4,i,j) = -CD(1)/abs(S(4,i,j));
    
    end
end


end
