% Jarikre Efe Jeffery
% 101008461
% Elec_4700 Assignment 1

clc
clearvars
close all
set(0,'DefaultFigureWindowStyle','docked')

m0 = 9.10938215e-31; % rest mass
effective_m = 0.26*m0; %effective mass
q0 = 1.60217653e-19;   %electron charge
dirac_c = 1.054571596e-34; %dirac's constant
planck_c = dirac_c*2*pi;     % Planck's constant
Boltz_c =  1.3806504e-23; % Boltzmann Constant
eps_0 = 8.854187817e-12; % permittivity of freespace
mu_0 = 1.2566370614e-6; % permeability of freespace
c = 299792458; % Speed of light
g = 9.80665; % gravityPa
am = 1.66053892e-27;

Thermal_v = sqrt((Boltz_c * 300) / effective_m); % Calculation for the thermal velocity
Stand_dev = Thermal_v/(sqrt(2)); % Standard Deviation of the two velocities

Time_bet_collisions = 0.2 * 10 ^ -12; % time between collisions
Time_step = 1000; % Time step
delta_t = 7.5 * 10 ^ -15; % Change in psoition of particles

wid_x = 200 * 10 ^ -9; % x-boundaries
len_y = 100 * 10 ^ -9; % y-boundaries.
size_particle = 100; % Number of particles
x_pos =  rand(1, size_particle) .* wid_x; % random x_position
y_pos =   rand(1, size_particle) .* len_y; % random y_position
theta = rand(1,size_particle)*2*pi;
x_vel = Thermal_v*cos(theta); % Assigning a random velocity to particles
y_vel = Thermal_v*sin(theta); % Assigning a random velocity to particles

figure (1)
xlabel("position of x");
ylabel("position of y");

isinbx = true;
while isinbx == true
    inbx = ((x_pos <= (1.15 * wid_x/2) & (x_pos >= (0.85 * wid_x/2))) & ((y_pos < (len_y/3)) | y_pos >= (2*len_y/3)));
    if (sum(inbx) > 0)
        x_pos(inbx) = rand(1, sum(inbx)) .* wid_x;
        y_pos(inbx) = rand(1, sum(inbx)) .* len_y;
    else
        isinbx = false;
        %break;
    end
   
end

% Implement a for loop here

for i = 1:200
    x_old = x_pos;
    x_pos = x_old + delta_t*x_vel;
    y_old = y_pos;
    y_pos = y_old + delta_t*y_vel;
   
    % Change in Velocity for the y-boundary condition
   
    index1 = (y_pos > len_y);
    index2 = (y_pos < 0);
    y_vel(index1) = -y_vel(index1);
    y_vel(index2) = -y_vel(index2);
   
   
    % implementing wrap around for the x-boundary condition
   
    index3 = ((x_pos) > wid_x); % Right ahnd position
    index4 = ((x_pos) < 0); % left hand position
   
    x_old(index3) = x_old(index3) - wid_x;
    x_pos(index3) = x_pos(index3) - wid_x;
   
    x_old(index4) = x_old(index4) + wid_x;
    x_pos(index4) = x_pos(index4) + wid_x;
   
    % plot particles moving
    plot (x_pos,y_pos,'.');
    pause(0.05)
    hold on
   
    % Draws the rectangular blocks
    line([0.85*wid_x/2 0.85*wid_x/2], [len_y 2*len_y/3]);
    line([1.15*wid_x/2 1.15*wid_x/2], [len_y 2*len_y/3]);
    line([0.85*wid_x/2 1.15*wid_x/2], [len_y len_y]);
    line([0.85*wid_x/2 1.15*wid_x/2], [2*len_y/3 2*len_y/3]);
   
    line([0.85*wid_x/2 0.85*wid_x/2], [0 len_y/3]);
    line([1.15*wid_x/2 1.15*wid_x/2], [0 len_y/3]);
    line([0.85*wid_x/2 1.15*wid_x/2], [0 0]);
    line([0.85*wid_x/2 1.15*wid_x/2], [len_y/3 len_y/3]);
    
    % Models regions where the particles bounce off the blocks
    inbx = ((x_pos <= (1.15 * wid_x/2) & (x_pos >= (0.85 * wid_x/2))) & ((y_pos < (len_y/3)) | y_pos >= (2*len_y/3)));
    between = ((x_old <= (1.15 * wid_x/2) & (x_old >= (0.85 * wid_x/2))) & ((y_old > (len_y/3)) & y_old <= (2*len_y/3)));
   
    % implements the change in velocities of the blocks when they hit the
    % rectangular blocks
    x_vel(inbx&(~between)) = -x_vel(inbx&(~between));
   
    y_vel(inbx&between) = -y_vel(inbx&between);
 
   
end