%% Pierce A Garver 1 November 2024
clear all
close all
%% set directories and import data
data_path = '';
filename = 'GloTEC_TEC_2024_11_21.nc';
nav_data_filename = 'full_states.csv';
vars = ncinfo([data_path filename]);
TEC = ncread([data_path filename], 'TEC'); % I think the TEC data is stored such that the x-dimension is latitude and y-dimension is longitude
long = ncread([data_path filename],'longitude');
lat = ncread([data_path filename], 'latitude');
collect_time = datetime('2024-11-21 12:57:35'); % Finds unix time for 18:20 GMT. This is the same as 10:20 PST.
collect_time.TimeZone = 'America/Los_Angeles';
collect_time_unix = posixtime(collect_time);
time = ncread([data_path filename], 'time'); %Time is reported in unix time
rounded_time = round_unix_time(time,collect_time_unix);
collect_time_index = find(rounded_time == time);
time_dmy = datetime(time(1,1),'ConvertFrom', 'posixtime'); % Convert unix time to GMT

%% Executing Newton Raphson to get rough initial estimate. NOTE: This code currently uses HW 2 Problem 1 Data.THIS WILL HAVE TO BE CHANGED
x_est = [0;0;0]; % Initial Estimate of position. Must be (3 x 1) vector
b_est = 0; % Clock bias estimate, should stay as 0 since psuedoranges are already corrected.
max_iterations = 100;
convergence_crit = .01; % Convergence is based on the 2-norm of the update
data_rows = 10176; % Number of rows in the given excel spreadsheet.
[pos_mat_array, corr_pr_vec_array] = read_data(strcat(data_path,nav_data_filename),data_rows);
alpha = [2.887*10^(-8) 0.0 -1.192*10^(-7) 1.788*10^(-7)];
beta = [1.434*10^5 -1.475*10^5 0 -1.311*10^5];
wgs84 = wgs84Ellipsoid('meter');
data_points = 173;%494;
for k = 1:data_points%length(pos_mat_array) 173
    k
    pos_mat = pos_mat_array{k}';
    corr_pr_vec = corr_pr_vec_array{k};
    for i =1:1
        %[pos_mat, corr_pr_vec] = read_data(strcat(data_path,nav_data_filename),data_rows(i));
        %pos_mat
        [opt_pos, opt_b, update, num_iter] = newton_raphson(x_est,b_est, pos_mat, corr_pr_vec, max_iterations, convergence_crit);
        opt_pos;
        opt_b;
        num_iter;
        update_mag = norm(update);
        opt_pos_mag = norm(opt_pos);
        user_lla = ecef2lla(opt_pos','WGS84');
    end
    for i =1:1
        corr_pr_vec = corr_pr_vec - opt_b;
        [opt_pos, opt_b, update, num_iter] = newton_raphson(x_est,b_est, pos_mat, corr_pr_vec, max_iterations, convergence_crit);
        opt_pos;
        opt_b;
        num_iter;
        update_mag = norm(update);
        opt_pos_mag = norm(opt_pos);
        user_lla = ecef2lla(opt_pos','WGS84');
    end

    for i =1:1
        dIon = [];
        for j = 1:size(pos_mat,2)
            sat_pos = pos_mat(:,i);
            [az, elev, slantRange] = ecef2aer(sat_pos(1),sat_pos(2),sat_pos(3),user_lla(1),user_lla(2),user_lla(3),wgs84);
            [week, tow] = cal2gpstime(2024, 11, 21, 12,57, 35);
            dIon(j) = klobuchar(user_lla(1),user_lla(2),elev,az,tow,alpha,beta);
        end
        corr_pr_vec_klob = corr_pr_vec - dIon';
        [opt_pos_klob, opt_b, update, num_iter] = newton_raphson(x_est,b_est, pos_mat, corr_pr_vec, max_iterations, convergence_crit);
        opt_pos;
        opt_b;
        klob_lla = ecef2lla(opt_pos_klob','WGS84');
    end
    
    %% Getting Lat and Lon of the Piercing Points
    %wgs84 = wgs84Ellipsoid('kilometer');
    pierce_point_lla_mat = [];
    for i = 1:size(pos_mat,2)
            sat_pos = pos_mat(:,i);
            [az, elev, slantRange] = ecef2aer(sat_pos(1),sat_pos(2),sat_pos(3),user_lla(1),user_lla(2),user_lla(3),wgs84);
            az_el_sR_mat(:,i) = [az; elev; slantRange]; % Creating a matrix which holds the az,el,SR vector from user to each satellite
            f2_layer_alt = 350000; %[m] F2 layer altitude. In reality it varies a lot but I don't think there is live data on F2 layer altitude 
            pierce_point_x_dist = f2_layer_alt/tand(elev);
            pierce_point_enu = [cosd(az)*pierce_point_x_dist; sind(az)*pierce_point_x_dist; f2_layer_alt];
            %pierce_point_enu = [sind(az)*cosd(elev)*pierce_point_x_dist; cosd(az)*cosd(elev)*pierce_point_x_dist; f2_layer_alt];
            pierce_point_lla = [user_lla(1)+pierce_point_enu(2)*(360/(40*10^6)); user_lla(2)+pierce_point_enu(1)*0.00001/cosd(user_lla(1)); f2_layer_alt];
            pierce_point_lla_mat(:,i) = pierce_point_lla;
    end
    
    %% Get TEC Values at the Pierceing Points for all of the Satellites in view
    receiver_lat = user_lla(1);
    receiver_lon = user_lla(2);
    for i = 1:size(pierce_point_lla_mat,2)
        [rounded_rec_lon, rounded_rec_lat] = round_coords(pierce_point_lla_mat(2,i), pierce_point_lla_mat(1,i));
        lat_index = find(lat==rounded_rec_lat);
        lon_index = find(long==rounded_rec_lon);
        TEC_value = TEC(lat_index,lon_index,collect_time_index);
        pierce_point_TEC_mat(i) = TEC_value;
    end
    
    %% Solving for ionosphere index of refraction
    h_iono_num = 500000; %[m] This is the thickness of the ionosphere. Very rough number which in fact changes all the time. Number determined from electron density plot fig 7.1 of "Introduction to the Space Environment second edition" by Thomas F. Tascione
    % Secondary note on h_iono: using the daytime night time curve, this h_iono
    % corresponds with the bottom being at 100km altitude where electron
    % density is about 10^5 [cm^-1] and the top being at 600 km  wehre electron
    ep_0_num = 8.854187817*10^(-12);
    m_e_num =  9.1093837*10^(-31); %[kg]
    q_num = 1.602176634*10^(-19); %[c]
    freq_num = 1575.42*10^6 * 2*pi; %[rad/sec] carrier frequency of the gps L1 signal
    syms n_ion TEC_sym h_iono ep_0 m_e q omega
    for i = 1:size(pierce_point_TEC_mat,2)
        omega_plasma = ((TEC_sym*q_num^2)/(ep_0_num*m_e_num))^(1/2);
        omega_plasma = double(subs(omega_plasma,TEC_sym, pierce_point_TEC_mat(i)*10^16/h_iono_num));
        N_ion = (1-omega_plasma^2/freq_num^2)^(1/2);
        iono_index_refraction(i) = N_ion; 
    end
    
    %% Solve for r1 and phi1
    syms r1 phi1
    h_tropo = 100000; %[m] height of the troposphere
    r1_vec = [];
    phi1_vec = [];
    for i = 1:length(corr_pr_vec)
        eq2 = r1*cosd(phi1) + h_iono_num + h_tropo == corr_pr_vec(i)*sind(az_el_sR_mat(2,i));
        eq3 = r1*sind(phi1) + h_iono_num*tand(asind(sind(phi1)/iono_index_refraction(i))) + h_tropo*tand(phi1) == corr_pr_vec(i)*cosd(az_el_sR_mat(2,i));
        [r1_num, phi1_num] = vpasolve([eq2 eq3], [r1 phi1]);
        r1_vec(i) = r1_num;
        phi1_vec(i) = phi1_num;
    end
    
    %% Calculate ionospheric correction time for each satellite
    c = 299792458;
    T_iono_corr_vec = [];
    iono_corr_vec_m = [];
    for i = 1:length(r1_vec)
        n_g = iono_index_refraction(i) + freq_num*double(subs(diff((1-(((TEC_sym*q_num^2)/(ep_0_num*m_e_num))^(1/2))^2/omega^2)^(1/2),omega), [TEC_sym omega], [pierce_point_TEC_mat(i)*10^16/h_iono_num freq_num]));
        T_iono_corr = double(((r1_vec(i)+sqrt(h_iono_num^2 + (h_iono_num*tand(asind(sind(phi1_vec(i)/iono_index_refraction(i)))))^2) + sqrt(h_tropo^2 + (h_tropo*tand(phi1_vec(i)))^2))-corr_pr_vec(i))/(c*n_g));
        T_iono_corr_vec(i) = T_iono_corr;
        iono_corr_vec_m(i) = T_iono_corr * c * n_g;
    end
    bad_entry_indices = find(abs(iono_corr_vec_m) > 500);
    iono_corr_vec_m(bad_entry_indices) = 0;

    %% Running Newton Raphson Again
    for i =1:1
        corr_pr_vec = corr_pr_vec - iono_corr_vec_m';
        [opt_pos_iono_corr, opt_b, update, num_iter] = newton_raphson(x_est,b_est, pos_mat, corr_pr_vec, max_iterations, convergence_crit);
        opt_pos;
        opt_b;
        num_iter;
        update_mag = norm(update);
        opt_pos_mag = norm(opt_pos_iono_corr);
        user_lla_iono_corrected = ecef2lla(opt_pos_iono_corr','WGS84');
    end
    ground_truth_ecef = lla2ecef([37.425556 -122.161389 0],'WGS84');
    iono_corr_error(k) = norm(opt_pos_iono_corr - ground_truth_ecef);
    uncor_error(k) = norm(opt_pos - ground_truth_ecef);
    klob_error(k) = norm(opt_pos_klob - ground_truth_ecef);
    user_lla_uncor_vec(:,k) = user_lla;
    user_lla_corr_vec(:,k) = user_lla_iono_corrected;
    user_lla_klob_vec(:,k) = klob_lla;
    
    
end
%%
klob_alt_error = user_lla_klob_vec(3,:) - 29;
iono_corr_alt_error = user_lla_corr_vec(3,:) - 29;
%% Plotting
fig1 = figure;
scatter(user_lla_uncor_vec(1,:),user_lla_uncor_vec(2,:),'filled');
hold on
scatter(user_lla_corr_vec(1,:),user_lla_corr_vec(2,:));
scatter(user_lla_klob_vec(1,:),user_lla_klob_vec(2,:),'green');
scatter(37.425556,-122.161389,'filled')
xlabel('Lattitude')
ylabel('Longitude')
title('Position Measurements')
legend('Uncorrected','Ionosphere Corrected','Klobuchar Corrected','Ground Truth')
grid on
hold off

fig2 = figure;
plotx = linspace(1,data_points,data_points);
plot(plotx,iono_corr_error)
hold on
plot(plotx,uncor_error)
plot(plotx,klob_error,'.-')
xlabel('Measurement')
ylabel('Error [m]')
legend('Corrected Data', 'Uncorrected Data', 'Klobuchar Model Error')
grid on
title('Errors')
hold off

fig3 = figure;  % This plot actually won't work yet because we don't have the ground truth altitude.
plot(plotx,klob_alt_error)
hold on
plot(plotx,iono_corr_alt_error)
grid on
xlabel('Measurement')
ylabel('Altitude Error [m]')
title('Altitude Errors')
legend('Klobuchar Corrected Alt Error','Iono Corrected Alt Error')
hold off

fig4 = figure;
plotx = linspace(1,data_points,data_points);
plot(plotx, klob_error - iono_corr_error)
hold on
xlabel('Measurement')
ylabel('Difference in Error [m]')
grid on
title('Difference of Errors [Klobuchar Corr Err - Homemade Corr Err')
hold off
%% FUCNTIONS


% round_coords rounds the receiver coordinates to the closest lat-long
% value stored in the lat-long variables of the .nc data file. This will be
% used for indexing so that we can select the correct TEC value for the
% receiver.
function [rounded_rec_lon, rounded_rec_lat] = round_coords(receiver_lon, receiver_lat)
    remainder = rem(abs(receiver_lon-2.5),5);
    lon_sign = (receiver_lon/abs(receiver_lon));
    lon_mult = floorDiv(abs(receiver_lon),5);
    if remainder > 5/2
        lon_mult = lon_mult+1;
    end
    rounded_rec_lon = (5*lon_mult+2.5)*lon_sign;
    remainder = rem(abs(receiver_lat)-1.25,2.5);
    lat_sign = (receiver_lat/abs(receiver_lat));
    lat_mult = floorDiv(abs(receiver_lat),2.5);
    if remainder > 2.5/2
        lat_mult = lat_mult+1;
    end
    rounded_rec_lat = (2.5*lat_mult+1.25)*lat_sign;
end

% Read and store data from an excel file for main problem.
% Loads the data and first gets rid of any emtpy rows. Then it chunks the
% data into time steps. If the GNSS logger app read the satellite signals
% at 200 different times, then there are 200 measurements. Each
% "measurement" includes the psuedorange to each satellite, and the
% satellite's position in ECEF. Each measurement is then stored as an
% array element inside of a larger cell array. The lenght of the cell array
% is the number of measurements.
function [pos_mat, corr_pr_vec] = read_data(filename, file_rows)
    data = readmatrix(filename, "Range", strcat('AN2:AX',string(file_rows)));
    empty_row_index = find(isnan(data(:,5)));
    data(empty_row_index,:) = [];
    i_start = 0;
    j=0;
    while i_start + 1 < size(data,1)-1
        j = j + 1;
        gps_millis = data(i_start+1,1);
        time_period_length = find(data(i_start+1:end,1) ~= gps_millis,1,'first')-1;
        pos_mat_page = data(i_start+1:time_period_length+i_start,5:7);
        pr_vec_page = data(i_start+1:time_period_length+i_start,2);
        sat_clock_bias = data(i_start+1:time_period_length+i_start,10);
        pr_vec_page = pr_vec_page + sat_clock_bias; % Correcting for satellite clock bias
        i_start = time_period_length + i_start;
        pos_mat{j} = pos_mat_page;
        corr_pr_vec{j} = pr_vec_page;
        i_start;
    end
    pos_mat(end)=[];
    corr_pr_vec(end)=[];
end

% Executes entire newton raphson process and iteration.
% Uses other helper functions to execute the newton raphson process and continues iteration until convergence is acheived or the maximum number of iterations is hit.
function [opt_pos, opt_b, update, num_iter] = newton_raphson(x_est, b_est, pos_mat, corr_pr_vec, max_iterations, convergence_crit) %Add opt clock bias as an output?
    num_iter = 0;
    update = convergence_crit+1;
    while norm(update) > convergence_crit && num_iter < max_iterations
        G = get_geometry_matrix(x_est, pos_mat);
        theory_prange_vec = get_theoretical_psuedoranges(x_est, b_est, pos_mat);
        delta_rho = theory_measured_dif(theory_prange_vec, corr_pr_vec);
        update = newton_raphson_step(G, delta_rho);
        [x_est, b_est] = update_current_pos(x_est, b_est, update);
        num_iter = num_iter + 1;
        %norm(update)
    end
    opt_pos = x_est;
    opt_b = b_est;
end

% Calculates Geometry Matrix
% Uses the position matrix for all satellites (3 x n) and the initial / current guess (3 x 1)  to calculate the Geometry matrix (n x 3), where n is the number of satellite signals received. Each row of the geometry matrix is a unit vector pointing in the direction of least error for that satellite.
function G = get_geometry_matrix(x_est, pos_mat)
    length = size(pos_mat,2);
    for i = 1:length
        pos_vec = pos_mat(:,i);
        unit_cor_vec = ((pos_vec - x_est)/norm(pos_vec - x_est))';
        G(i,:) = unit_cor_vec;
    end
    G(:,4) = ones(length,1);
end

% Calculate theoretical psuedorange measurement
% Takes the current estimate of position (3 x 1) and clock bias (scalar) along with the pos_matrix of satellites. It then calculates the difference between each satellite's position vector and the estimate and calculates the Euclidean norm. Then, all of these 2-norms are stored in a (n x 1) vector.
function theory_prange_vec = get_theoretical_psuedoranges(x_est, b_est, pos_mat)
    length = size(pos_mat,2);
    for i = 1:length
        pos_vec = pos_mat(:,i);
        theory_prange = norm(pos_vec-x_est) + b_est;
        theory_prange_vec(i,1) = theory_prange;
    end
end

% Calculates difference between theoretical pranges and measured pranges. Returns delta_rho which is a (n x 1) vector.
function delta_rho = theory_measured_dif(theory, measured)
 delta_rho = theory - measured;
end

% Calculate update for one step of newton raphson
% Executes the matrix math to calculate the update (3 x 1) using the geometry matrix and delta_rho vector.
function update = newton_raphson_step(G, delta_rho)
    update = inv(G'*G)*G'*delta_rho;
end

% Calculate the new position estimate.
% Adds the update to the current position to move one step closer to the least squares solution.
function [new_est, new_b] = update_current_pos(current_pos, current_b, update)
    new_est = current_pos + update(1:3);
    new_b = current_b - update(4);
end

function rounded_time = round_unix_time(unix_time_vec, unrounded_unix_time)
    rounded_time = interp1(double(unix_time_vec), double(unix_time_vec), unrounded_unix_time, 'nearest');
end