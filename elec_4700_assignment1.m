% Jarikre Efe Jeffery
% 101008461
% Elec_4700 Assignment 1

clc
clearvars
close all
set(0,'DefaultFigureWindowStyle','docked')

simulation_time = 100;
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
temperature = 300; % Room temperature value
am = 1.66053892e-27;

Thermal_v = sqrt((Boltz_c * 300) / effective_m); % Calculation for the thermal velocity
stand_dev = Thermal_v/(sqrt(2)); % Standard Deviation of the two velocities

time_bet_collisions = 0.2 * 10 ^ -12; % time between collisions
time_step = 1e-15; % Time step
delta_t = 7.5 * 10 ^ -15; % Change in psoition of particles
wid_x = 200 * 10 ^ -9; % x-boundaries
len_y = 100 * 10 ^ -9; % y-boundaries.
size_particle = 50; % Number of particles
init_xv = randn(1,size_particle) .*stand_dev;
init_yv = randn(1,size_particle) .*stand_dev;
x_pos =  rand(1, size_particle) .* wid_x; % random x_position
y_pos =   rand(1, size_particle) .* len_y; % random y_position
theta = rand(1,size_particle)*2*pi;

sigma = sqrt(Boltz_c*temperature/effective_m)/4;
mu = Thermal_v;
MB_dist = makedist('Normal',mu,sigma);
vel = random(MB_dist,1,size_particle);

x_vel = Thermal_v*cos(theta); % Assigning a random velocity to particles
y_vel = Thermal_v*sin(theta); % Assigning a random velocity to particles

figure (1)
xlabel("position of x");
ylabel("position of y");

temparr = zeros(1, simulation_time);
tmpx = (1:1:simulation_time);

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

v_rms = sqrt((x_vel .^ 2) + (y_vel .^ 2));
v_rms_arr = zeros(1, 1000);

% Calculating the probability of scattering
pscat = 1 - (exp((-1 * delta_t) / (time_bet_collisions)));
index_is2 = zeros(1, size_particle);

numcol = 0;
temp_diff = 0;
sum_temp_diff = 0;

% Implement a for loop here
for i = 1:simulation_time
   
    index_is = pscat > rand(1,size_particle);
   
    if index_is(1) ~= index_is2(1)
        numcol = numcol + 1;
        sumtdfiff = sum_temp_diff + temp_diff;
        temp_diff = 0;
    else
        temp_diff = temp_diff + 1;
        %**************Use time time travelled times speed to find mfp
    end
   
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
   
    line([0.85*wid_x/2 0.85*wid_x/2], [len_y 2*len_y/3]);
    line([1.15*wid_x/2 1.15*wid_x/2], [len_y 2*len_y/3]);
    line([0.85*wid_x/2 1.15*wid_x/2], [len_y len_y]);
    line([0.85*wid_x/2 1.15*wid_x/2], [2*len_y/3 2*len_y/3]);
   
    line([0.85*wid_x/2 0.85*wid_x/2], [0 len_y/3]);
    line([1.15*wid_x/2 1.15*wid_x/2], [0 len_y/3]);
    line([0.85*wid_x/2 1.15*wid_x/2], [0 0]);
    line([0.85*wid_x/2 1.15*wid_x/2], [len_y/3 len_y/3]);
   
    inbx = ((x_pos <= (1.15 * wid_x/2) & (x_pos >= (0.85 * wid_x/2))) & ((y_pos < (len_y/3)) | y_pos >= (2*len_y/3)));
    between = ((x_old <= (1.15 * wid_x/2) & (x_old >= (0.85 * wid_x/2))) & ((y_old > (len_y/3)) & y_old <= (2*len_y/3)));
   
    x_vel(inbx&(~between)) = -x_vel(inbx&(~between));
   
    y_vel(inbx&between) = -y_vel(inbx&between);
   
    v_rms = sqrt ((init_xv.^2)+(init_yv.^2));
    show1 = (sqrt (2) * (mean(v_rms)^2) * effective_m )/Boltz_c;
    temparr (1,i)=  show1;
   
end

[x_graph, y_graph] = meshgrid(0:(wid_x/10):wid_x, 0:(len_y/10):len_y);
electron_mate = zeros(11, 11);
temperature_emate = zeros(11, 11);
numelectrol = 0;
totalvel = 0;

for j = 1:10
    x_min = x_graph(1, j);
    x_max = x_graph(1,j+1);
    for k = 1:10
        y_min = y_graph(k, 1);
        y_max = y_graph(k+1, 1);
        for m = 1:size_particle
            if((x_pos(m) > x_min) && (x_pos(m) < x_max) && ((y_pos(m) > y_min) && y_pos(m) < y_max))
                numelectrol = numelectrol + 1;
                electron_mate(j, k) = electron_mate(j, j) + 1;
                totalvel = totalvel + sqrt((x_vel(m) .^ 2) + (x_vel(m) .^ 2));
                temperature_emate(j, k) = ((sqrt(2)*(totalvel/numelectrol) ^ 2) * effective_m) / Boltz_c;
            end
        end
        totalvel = 0;
        numelectrol = 0;
    end
end

figure (2)
plot(tmpx,temparr)
title(["Graph of Temperature against time: Avg temp = " num2str(show1)]);
xlabel("time");
ylabel("Temperature");
hold on

% Mean Free Path and Mean Free Time Calculations
mean_free_time = (sum_temp_diff * delta_t)/numcol;
mean_free_path = (mean(v_rms) * mean_free_time);

fprintf("The Mean Free Time is = %12.15f\n", mean_free_time);
fprintf("The Mean Free Path is = %12.15f\n", mean_free_path);

figure(3); histogram(vel); title('Histogram of Thermal Velocities');
figure(4); surf(electron_mate); title('Density Mapping');
figure(5); surf(temperature_emate); title('Temperature Mapping');
end
