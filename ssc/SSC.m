%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSC.m - WFSim-based Controller for SOWFA
%
% This script implements a controller that:
% 1. Initializes WFSim for the 9-turbine layout
% 2. Receives measurements from SOWFA via ZeroMQ
% 3. Uses WFSim to compute flow field and sensitivities
% 4. Uses SOWFA power measurements for gradient-based optimization
% 5. Converts axial induction factors to pitch/torque/yaw commands
% 6. Sends control commands back to SOWFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;
disp('Starting WFSim-based wind farm controller for SOWFA...');


%% 1. Setup zeroMQ server
zmqServer = zeromqObj('/home/shijiehuang/OpenFOAM/zeroMQ/jeromq/target/jeromq-0.4.3.jar',1729,7200,true);
disp('ZeroMQ server started on port 4120');

%% 2. Initialize WFSim
% Add WFSim to path
addpath(genpath('/home/shijiehuang/wfcontrol_sj/WFSim'));

% Initialize the 9-turbine laayou in WFSim
Wp = layoutSet_sowfa_9turb_apc_alm_turbl();

% Set WFSim opions
modelOptions.Linearversion = 1;    % Linearized WFSim version
modelOptions.Projection = 0;       % No projection
modelOptions.exportLinearSol = 0;  % Don't export linear solution
modelOptions.exportPressures = 0;  % Don't export pressures
modelOptions.printConvergence = 0; % Don't print convergence
modelOptions.conv_eps = 1e-6;      % Convergence threshold
modelOptions.max_it_dyn = 1;       % Max iterations

% Initialize WFSim
[Wp, sol, sys] = InitWFSim(Wp, modelOptions, 0);
disp('WFSim initialized successfully');

%% 3. Initialize control variables
nTurbs = Wp.turbine.N;

% Initial control settings
CT_prime = 2 * ones(nTurbs, 1); % starting from greedy control might cause singular behaviour
yaw_angle = zeros(nTurbs, 1);

% Control constraints
CT_prime_min = 0.4 * ones(nTurbs, 1);
CT_prime_max = 3.6 * ones(nTurbs, 1);
yaw_min = -30 * ones(nTurbs, 1);
yaw_max	= 30 *	ones(nTurbs, 1);

% Controller parameters
alpha = 0.3;			% Learning rate for CT_prime (conservative for stability) 
alpha_Phi = 3; 		% Learning rate for yaw
mu = 0.00028;			% Regularization for CT_prime
mu_Phi = 0.0001; 		% Regularization for yaw
power_ref = nTurbs * 2e6;

% For conversion to SOWFA inputs
Kopt = 0.0255764; 		% Default NREL 5MW k-optimal
aOptimal = 1/3;			% Betz-optimal axial induction factor

%% 4. Controller main loop
counter = 0;
disp('Entering wind farm controller loop...');

% Store power and control
powerHistory = [];
CTHistory = [];
yawHistory = [];
windSpeedHistory = [];

% SOWFA time information
currentTime = 0;
lastOptimizationTime = 0;
optimizationInterval = 45; 	% Optimize every 60 seconds (Without such an interval, you'd be optimizing based on measurements that don't yet reflect your last control action) 

% Convergence tracking variables
powerImprovements = [];
lastPowerTotal = 0;

% Write initial settings to log file
logFile = fopen('controller_log.txt', 'w');
fprintf(logFile, 'WFSim-Based Controller for 9-Turbine SOWFA Simulation (ADM model)\n');
fprintf(logFile, 'Time: 0, Initial CT_prime: %f, Initial yaw: %f\n', CT_prime(1),yaw_angle(1));
fclose(logFile);

disp('IMPORTANT: Using relationship calibrated with ALMAdvanced model');
disp('K/K_opt = 24.6850*a^2 - 21.3041*a + 5.1988');
disp('Note that differences may exist when applied to ADM model');

% Control mode (0: greedy baseline, 1: WFSim-based)
% controlMode = 0;		% Start with greedy control as baseline
% switchTime = 300; 		% Switch to WFSim-based control after 300 seconds
% disp(['Starting with greedy control, switching to WFSim-based at ' num2str(switchTime) ' seconds']); 
disp('Starting directly with WFSim-based control strategy');
while 1
    % 1. Receive information from SOWFA
    try
        dataReceived = zmqServer.receive();
    catch e
        disp(['Error receiving data: ' e.message]);
        break;
    end

    currentTime  = dataReceived(1,1);
    measurementVector = dataReceived(1,2:end); % [powerGenerator[1], torqueRotor[1], thrust[1], powerGenerator[2], torqueRotor[2], thrust[2]]
    
    % 2. Extract measurements: [genPower,rotSpeedF,azimuth,rotThrust,rotTorque,genTorque,nacYaw,bladePitch]
    generatorPowerArray = measurementVector(1:8:end);
    rotorSpeedArray     = measurementVector(2:8:end);
    azimuthAngleArray   = measurementVector(3:8:end);
    rotorThrustArray    = measurementVector(4:8:end);
    rotorTorqueArray    = measurementVector(5:8:end);
    genTorqueArray      = measurementVector(6:8:end);
    nacelleYawArray     = measurementVector(7:8:end);
    bladePitchArray     = measurementVector(8:8:end);
    
    % Store power data from SOWFA measurements
    powerHistory = [powerHistory; generatorPowerArray'];

    % Log the current measurements periodically (every 10 seconds to avoid excessive output)
    if mod(counter, 10) == 0
	totalPower = sum(generatorPowerArray);
        disp(['Time: ' num2str(currentTime) ' s | Total power: ' num2str(totalPower/1e6) ' MW']);
        
        % Log individual turbine data
        for i = 1:nTurbs
            disp(sprintf('Turbine %d: Power = %.2f MW, RPM = %.1f', ...
                 i-1, generatorPowerArray(i)/1e6, rotorSpeedArray(i)));
        end

	% Append to log file
	logFile = fopen('controller_log.txt','a');
	fprintf(logFile, 'Time: %f, Total Power: %f MW\n', currentTime, totalPower/1e6);
	fclose(logFile);
    end

    % 3. Run optimization if enough time has passed
    if (currentTime - lastOptimizationTime >= optimizationInterval)
	lastOptimizationTime = currentTime;
	cycleCounter = floor(currentTime/optimizationInterval);
        disp(['Running feedback optimization cycle #' num2str(cycleCounter)]); 
    
	% 3.1 Use WFSim to simulate current state and linearize the model
	sol.k = sol.k + 1;
	[sol, sys] = Make_Ax_b_hsj(Wp, sys, sol, yaw_angle, CT_prime, modelOptions);
	[sol, sys] = Computesol(sys, sol, sol.k, modelOptions);
	[sol, eps] = MapSolution(Wp, sol, sol.k, modelOptions);

	% 3.2 Extract flow field and turbine inflow speeds
	uss = sol.u;
	linear_model = struct;
	linear_model.ss.U = zeros(nTurbs, 1);

	for j=1:nTurbs
	    linear_model.ss.U(j) = sum(uss(Wp.mesh.xline(j),Wp.mesh.yline{j}))/length(Wp.mesh.yline{j});
	end

	% Store estimated wind speeds
	windSpeedHistory = [windSpeedHistory; linear_model.ss.U'];
	
	% 3.3 Linearize the model
        try
            linear_model.A = sys.A\sys.Al;
            linear_model.B = sys.A\sys.Bl;
            linear_model.BPhi = sys.A\sys.BlPhi;

            nx = size(linear_model.A,1);

            % 3.4 Setup measurement equations
            linear_model.Cz = zeros(nTurbs, size(sys.A,1));
            for j=1:nTurbs
                linear_model.Cz(j,(Wp.mesh.xline(j)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{j}(1)-1:...
                    (Wp.mesh.xline(j)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{j}(end)-1) = 1/length(Wp.mesh.yline{j});
            end

            % 3.5 Setup power sensitivity matrices
            Rho = Wp.site.Rho;
            Drotor = Wp.turbine.Drotor;
            powerscale = Wp.turbine.powerscale;
            Ar = pi*(0.5*Drotor)^2;
            constant = powerscale*.5*Rho*Ar;

            linear_model.Cx = diag(3*constant*linear_model.ss.U.^2.*CT_prime.*(cos(yaw_angle*pi/180).^3));
            linear_model.C = linear_model.Cx * linear_model.Cz;
            linear_model.D = diag(constant*(linear_model.ss.U .* cos(yaw_angle*pi/180)).^3);
            linear_model.DPhi = diag(-1/180*pi*3*constant*linear_model.ss.U.^3.*CT_prime.*(cos(yaw_angle*pi/180)).^2.*sin(yaw_angle*pi/180));

            % 3.6 Steady-state mapping
            steady_matrix = linear_model.C * ((eye(nx) - linear_model.A)\linear_model.B) + linear_model.D;
            steady_matrix_Phi = linear_model.C * ((eye(nx) - linear_model.A)\linear_model.BPhi) + linear_model.DPhi;
            
            % 3.7 Compute gradients using SOWFA power measurements
            totalPower = sum(generatorPowerArray);  % Use SOWFA measurements
            nabla_CT_prime = 2 * ((steady_matrix)'*ones(nTurbs,1)*(totalPower - power_ref)/power_ref^2 + mu * CT_prime);
            nabla_yaw_angle = 2 * ((steady_matrix_Phi)'*ones(nTurbs,1)*(totalPower - power_ref)/power_ref^2 + mu_Phi * yaw_angle);
            
            % 3.8 Update control variables with projected gradient descent
            CT_prime_new = min([CT_prime_max, max([CT_prime_min, CT_prime - alpha * nabla_CT_prime],[],2)],[],2);
            yaw_angle_new = min([yaw_max, max([yaw_min, yaw_angle - alpha_Phi * nabla_yaw_angle],[],2)],[],2);
            
            % 3.9 Apply changes only if significant
            CT_prime_change = norm(CT_prime_new - CT_prime) / norm(CT_prime);
            yaw_angle_change = norm(yaw_angle_new - yaw_angle) / (norm(yaw_angle) + 1e-6);

            if CT_prime_change > 0.001 || yaw_angle_change > 0.001
                disp('Applying significant control updates');
                CT_prime = CT_prime_new;
                yaw_angle = yaw_angle_new;
                
                % Log control changes
                logFile = fopen('controller_log.txt', 'a');
                fprintf(logFile, 'Time: %f, Updated Controls - CT_prime=%f, Yaw=%f\n', currentTime, mean(CT_prime), mean(yaw_angle));
                fclose(logFile);
            else
                disp('Control updates below threshold, maintaining current values');
            end
        catch e
            disp(['Optimization error: ' e.message]);
            disp('Maintaining current control settings');
        end

        % Add convergence tracking
        totalPower = sum(generatorPowerArray);
        if lastPowerTotal > 0
            improvement = (totalPower - lastPowerTotal)/lastPowerTotal * 100;
            powerImprovements = [powerImprovements; improvement];
            
            % Enhanced logging
            disp(['Power improvement: ' num2str(improvement) '%']);
            
            logFile = fopen('controller_log.txt', 'a');
            fprintf(logFile, 'Cycle %d, Time: %f, Power: %f, Improvement: %f%%\n', ...
                cycleCounter, currentTime, totalPower/1e6, improvement);
            fclose(logFile);
        end
        lastPowerTotal = totalPower;
    end

    % 4. Store control history
    CTHistory = [CTHistory; CT_prime'];
    yawHistory = [yawHistory; yaw_angle'];
    
    % 5. Convert CT_prime to blade pitch and generator torque
    % Initialize control outputs
    torqueArrayOut = zeros(1, nTurbs);
    pitchAngleArrayOut = zeros(1, nTurbs);
    yawAngleArrayOut = 270 + yaw_angle'; % SOWFA uses 270 as aligned with flow

    for i = 1:nTurbs
        % 5.1 Convert CT_prime to axial induction factor
        % Using CT_prime = 4a/(1-a)
        % Solving for a: a = CT_prime/(4 + CT_prime)
        a = CT_prime(i)/(4 + CT_prime(i));
        
        % 5.2 All turbines are in Region 2 (8 m/s is below rated)
        % Set pitch to zero for maximum aerodynamic efficiency
        pitchAngleArrayOut(i) = 0.0;
        
        % 5.3 Get the k-scaling factor from your calibrated formula
        % K/K_opt = 24.6850*a^2 - 21.3041*a + 5.1988
        kScaling = 24.6850*a^2 - 21.3041*a + 5.1988;
        
        % Constrain scaling within reasonable limits
        kScaling = max(0.3, min(3.0, kScaling));
        
        % 5.4 Apply torque = k*omega^2 control law with scaled k
        rotorSpeedRadSec = rotorSpeedArray(i) * 2 * pi /60; % Convert rpm to rad/s
        kFactor = Kopt * kScaling;
        torqueArrayOut(i) = kFactor * rotorSpeedRadSec^2;
        
        % 5.5 Log for debug (only for first turbine, occasionally)
        if mod(counter, 100) == 0 && i == 1
            disp(['Turbine ' num2str(i) ': CT_prime=' num2str(CT_prime(i)) ...
                  ', a=' num2str(a) ...
                  ', kScaling=' num2str(kScaling) ...
                  ', torque=' num2str(torqueArrayOut(i)/1e3) ' kNm']);
        end
    end

    % 6. Create updated control signal string for SOWFA
    dataSend = setupZmqSignal(torqueArrayOut, yawAngleArrayOut, pitchAngleArrayOut);
    
    % 7. Send control actions back to SOWFA
    try
        zmqServer.send(dataSend);
    catch e
        disp(['Error sending data: ' e.message]);
        break;
    end
    
    % 8. Save control/power data periodically
    if mod(counter, 50) == 0
        try
            save('WFSim_control_data.mat', 'powerHistory', 'CTHistory', 'yawHistory', 'windSpeedHistory', 'currentTime', 'powerImprovements');
        catch
            disp('Warning: Unable to save control data');
        end
    end
    
    counter = counter + 1;
end

% Close connection and save final data
try
    zmqServer.disconnect();
    save('WFSim_control_data_final.mat', 'powerHistory', 'CTHistory', 'yawHistory', 'windSpeedHistory');
catch
    disp('Warning: Error during cleanup');
end

%% Helper Functions
% Function to create the control signal string for SOWFA
function [dataOut] = setupZmqSignal(torqueSignals, yawAngles, pitchAngles)
    dataOut = [];
    for i = 1:length(yawAngles)
        dataOut = [dataOut torqueSignals(i) yawAngles(i) pitchAngles(i)];
    end
end
