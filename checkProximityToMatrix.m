function checkProximityToMatrix(h0,v0,thetaPol1,thetaPol2,targetMat)
%% Plot the proximity of every matrix to the target matrix as a function
%% of HWP angles
%
% CHECKPROXIMITYTOMATRIX show the proximity of density matrices to the
% target matrix. The system consists of 2 sequential configurations of a
% calcite crystal, polarizer, HWP, where each of the HWP and polarizers can
% have a varying angle.
% 
%   CHECKPROXIMITYTOMATRIX(H,V,THETA1,THETA2, TARGETMAT) - Shows the proximity of the
%   density matrix for the input polarization going through the system
%   while the first HWP is in the angle represented by the x axis, and the
%   second HWP is in the angle represented by the y axis. A 1 value of a
%   given cell represents a 100% match between the density matrix and the
%   target matrix. The polarizers are set to theta1 and theta2 degrees
%   respectively (H and V transmitting).
%
%   CHECKPROXIMITYTOMATRIX(H,V,TARGETMAT) - Shows the proximity of the
%   density matrix for the input polarization going through the system
%   while the first HWP is in the angle represented by the x axis, and the
%   second HWP is in the angle represented by the y axis. A 1 value of a
%   given cell represents a 100% match between the density matrix and the
%   target matrix. The polarizers are automatically set to 0 and 90 degrees
%   respectively (H and V transmitting).
%
%   CHECKPROXIMITYTOMATRIX(H,V) - Shows the proximity of the
%   density matrix for the input polarization going through the system
%   while the first HWP is in the angle represented by the x axis, and the
%   second HWP is in the angle represented by the y axis. A 1 value of a
%   given cell represents a 100% match between the density matrix and the
%   target matrix. In this case the default choice for the target matrix is
%   an even matrix where every cell in the 6x6 matrix equals 1/6 exactly.
%   The polarizers are automatically set to 0 and 90 degrees respectively
%   (H and V transmitting).

%% Set default target Matrix
    % set default num of steps if not entered
    if ~exist('thetaPol1','var')
        pol1Angle = deg2rad(0);
    else
        pol1Angle = deg2rad(thetaPol1);
    end
    % set default polarizer 1 angle if not entered
    if ~exist('thetaPol2','var')
        pol2Angle = deg2rad(90);
    else
        pol2Angle = deg2rad(thetaPol2);
    end
    % set default polarizer 2 angle if not entered
    if ~exist('targetMat','var')
        targetMat = 1/6*ones(6);
    end
%% Create Data Matrix
    proximity = zeros(90);
    for firstTheta = 1:91
        for secTheta = 1:91
            currMat = calcDensityMatrix(h0,v0,firstTheta-1,secTheta-1,pol1Angle,pol2Angle);
            proximity(firstTheta,secTheta) = showFidelity(targetMat,currMat);
        end
    end
%% Display Heat Map of proximity
figure('Name','Plot of theta1 = 0 and theta2 = 1:90');
plot(0:90,proximity(1,:));
figure('Name','Heat map of angles proximity');
surf(0:90,0:90,proximity);
xlabel('\theta_1 - Angle of first HWP');
ylabel('\theta_2 - Angle of second HWP');
end

function fidelity = showFidelity(targetMat,compareMat)
%% Show Fidelity
% 
% SHOWFIDELITY shows how close the compareMat is to the targetMat, by
% showing the trace of their scalar multiplication.
%
%   SHOWFIDELITY(TARGETMAT,COMPAREMAT) - compare the compare matrix to the
%   target matrix.
    if (size(targetMat) ~= size(compareMat))
        error('Illegal input, matrix sizes must be same');
    end
    fidelity = trace(targetMat*compareMat);
end

function density = calcDensityMatrix(h0,v0,thetaHWP1,thetaHWP2,thetaPol1,thetaPol2)
    %% Check input validity
    % precision of input
    epsilon = 5e-16;
    % verify normalized input polarization
    assert(abs(abs(h0)^2+abs(v0)^2-1) < epsilon,'Input polarization is not normalized');
    % verify positive angle between 0 and 360
    assert((thetaHWP1 >= 0) && (thetaHWP1 <= 360),'Illegal angle, please input angle between 0 and 360');
    assert((thetaHWP2 >= 0) && (thetaHWP2 <= 360),'Illegal angle, please input angle between 0 and 360');
    assert((thetaPol1 >= 0) && (thetaPol1 <= 360),'Illegal angle, please input angle between 0 and 360');
    assert((thetaPol2 >= 0) && (thetaPol2 <= 360),'Illegal angle, please input angle between 0 and 360');
    
    %% Define Polarization and maximum angle to display
    % Polarity definition
    alpha = h0;
    beta = v0;

    %% Constants Definitions
    % Calcite Operator Matrix:
    U_cal = createCalMat();
    U_cal_perp = createPerpCalMat();
    % HWP Operator Matrix:
    U_hwp = createHWPMatrix(thetaHWP1);
    U_hwp_last = createHWPMatrix(-thetaHWP2);
    % HWP Operator Matrix:
    U_pol = createPolarMat(thetaPol1);
    U_pol_last = createPolarMat(thetaPol2);
    % Polarization Vector:
    polarization = zeros(6,1);
    polarization(1,1) = alpha;
    polarization(4,1) = beta;

    % Define the matrices of the HWPs per this angle, the output
    % polarization vector, and the density matrix.
    outputPolarization = U_hwp_last*U_pol_last*U_cal_perp*U_hwp*U_pol*U_cal*polarization;
    density = outputPolarization*(outputPolarization)';
end

% Create the half wave plate matrix representation
function U_hwp_mat = createHWPMatrix(theta)
    U_hwp_mat = zeros(2);
    U_hwp_mat(1,1) = cos(2*theta);
    U_hwp_mat(1,2) = sin(2*theta);
    U_hwp_mat(2,1) = sin(2*theta);
    U_hwp_mat(2,2) = -cos(2*theta);
    U_hwp_mat = kron(U_hwp_mat,eye(3));
end

% Create Polarizer Matrix in the angle theta
function U_polar_mat = createPolarMat(theta)
    U_polar_mat = zeros(2);
    U_polar_mat(1,1) = (cos(theta))^2;
    U_polar_mat(1,2) = sin(theta)*cos(theta);
    U_polar_mat(2,1) = sin(theta)*cos(theta);
    U_polar_mat(2,2) = (sin(theta))^2;
    U_polar_mat = kron(U_polar_mat,eye(3));
end

% Create Calcite Operator Matrix Representation Delays V
% polarization
function U_cal_mat = createCalMat()
    U_cal_mat = zeros(6);
    U_cal_mat(1,1) = 1;
    U_cal_mat(2,2) = 1;
    U_cal_mat(3,3) = 1;
    U_cal_mat(5,4) = 1;
    U_cal_mat(6,5) = 1;
end

% Create Perpendicular Calcite Operator Matrix Representation Delays H
% polarization
function U_perp_cal_mat = createPerpCalMat()
    U_perp_cal_mat = zeros(6);
    U_perp_cal_mat(4,4) = 1;
    U_perp_cal_mat(5,5) = 1;
    U_perp_cal_mat(6,6) = 1;
    U_perp_cal_mat(2,1) = 1;
    U_perp_cal_mat(3,2) = 1;
end