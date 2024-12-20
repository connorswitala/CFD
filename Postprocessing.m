close all;

Primitives = zeros(5,Nx,Ny);

%% Values for non-dimensionalizing viscous plate. 
Mach = 1.5;
Temp = 300;
gamma = 1.4;
R = 287;
sos = sqrt(gamma*R*Temp);
density = 0.01;
uvel = Mach*sos;
vvel = 0;
Pressure = density*R*Temp;

for h = 1:Nx
    for b = 1:Ny
       [Primitives(:,h,b)] = convert_to_primitives(U(:,h,b));
    end
end

xcenter = squeeze(center(1,:,:));
ycenter = squeeze(center(2,:,:));

%% Variables to plot

rho = squeeze(Primitives(1,:,:));
u = squeeze(Primitives(2,:,:));
v = squeeze(Primitives(3,:,:));
p = squeeze(Primitives(4,:,:));
T = squeeze(Primitives(4,:,:));
M = u./(sqrt(gamma.*R.*T));

%% For non-dimensionalized viscous plate plot
uint = u./uvel;
tint = T./Temp;


figure(1);
a = contourf(xcenter, ycenter, u, 'LineColor', 'none');
xlabel('x');
ylabel('y');
grid on;
grid minor;
colormap jet;
colorbar;
axis equal
% Pfilename = sprintf('Pressure_%d_by_%d_Ramp.png', Nx, Ny);
% print(gcf, Pfilename, '-dpng', '-r800');

% figure(2);
% a = semilogy(t,outer_residual,'-b');
% a.LineWidth = 1;
% grid on;
% grid minor;
% ylabel('R_i');
% xlabel('time (s)');
% Rfilename = sprintf('Residual_%d_by_%d_Ramp.png', Nx, Ny);
% print(gcf, Rfilename, '-dpng', '-r800');

% figure(3);
% a = contourf(xcenter,ycenter, M, 'LineColor', 'none');
% xlabel('x');
% ylabel('y');
% grid on;
% grid minor;
% colormap jet;
% colorbar;
% Mfilename = sprintf('Mach_%d_by_%d_Ramp.png', Nx, Ny);
% print(gcf, Mfilename, '-dpng', '-r800');

% figure(4);
% a = contourf(xcenter,ycenter,T, 'LineColor', 'none');
% xlabel('x');
% ylabel('y');
% grid on;
% grid minor;
% colormap jet;
% colorbar;
% % Tfilename = sprintf('Temperature_%d_by_%d_Ramp.png', Nx, Ny);
% % print(gcf, Tfilename, '-dpng', '-r800');

% figure(5);
% clf;
% a = plot(uint(Nx/2,:), ycenter(Nx/2,:), '-.r');
% a.LineWidth = 2;
% grid on;
% hold on;
% b = plot(tint(Nx/2,:), ycenter(Nx/2,:), '--b');
% b.LineWidth = 2;
% legend('Velocity', 'Temperature');
% xlabel('Non-dimensionalized ratio of temperature/velocity to free stream values.')
% ylabel('y');

%% Function to convert U to primitive variables

function [W] = convert_to_primitives(U)
           
    gamma = 1.4;    
    R = 287;
    % Extract conserved variables from U
    rho = U(1);          % Density
    rhou = U(2);         % rho * u
    rhov = U(3);         % rho * v
    E = U(4);            % Total energy

    % Calculate primitive variables
    W = zeros(5,1);
    W(1) = rho;                  % Density is directly available
    W(2) = rhou / rho;          % Velocity in x-direction: u = (rho * u) / rho
    W(3) = rhov / rho;          % Velocity in y-direction: v = (rho * v) / rho

    % Total energy equation to calculate pressure
    kinetic_energy = 0.5 * rho * (W(2)^2 + W(3)^2);
    W(4) = (gamma - 1) * (E - kinetic_energy);  % Pressure: p = (gamma - 1) * (E - 0.5*rho*(u^2 + v^2))
    W(5) = W(4)/(R*W(1));

end
