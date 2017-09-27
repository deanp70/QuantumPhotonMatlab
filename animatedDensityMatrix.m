function dens = animatedDensityMatrix(h0,v0,theta,numOfSteps)
%% Animated Density Matrix of PPT (3 degrees of freedom in time axis)
%
% ANIMATEDDENSITYMATRIX animate a density matrix in a system of two 
% calcite crystals with a half wave plate of variable theta between 
% them and another half wave plate of minus theta after the second crystal.
%
%   ANIMATEDDENSITYMATRIX(H,V,THETA) - displays the density matrix for
%   the HWP positioned at theta (degrees) for the input polarization 
%   (H,V).
%
%   ANIMATEDDENSITYMATRIX(H,V,THETA,STEPNUM) - simulates the turning
%   of the HWP between 0 and theta (degrees) for the input polarization
%   (H,V) in STEPNUM steps (recommend using more than 100 steps at least
%   for fluid animation). 
%   Keyboard based ui:
%   press 'q' - to exit.
%   press 'o' - to move back to the previous angle.
%   press 'p' - to move forward to the next angle.
%

    %% Check input validity
    % precision of input
    epsilon = 5e-16;

    % validate number of args
    if (nargin < 3) || (nargin > 5)
        error('Illegal number of arguments, 3 or 4 args expected');
    end
    % verify normalized input polarization
    assert(abs(abs(h0)^2+abs(v0)^2-1) < epsilon,'Input polarization is not normalized');
    % verify positive angle between 0 and 360
    assert((theta >= 0) && (theta <= 360),'Illegal angle, please input angle between 0 and 360');
    
    %% Define Polarization and maximum angle to display
    % Polarity definition
    alpha = h0;
    beta = v0;

    % Final Angle of the simulation
    endAngle = deg2rad(theta);

    %% Constants Definitions
    % First 0 angle HWP
    % U_hwp_first = createHWPMatrix(0);
    % Last 22.5 angle HWP
    % U_hwp_last = createHWPMatrix(pi/4);
    % Create Polarization vector (in 6 dimensional space)
    polarization = zeros(6,1);
    polarization(1,1) = alpha;
    polarization(4,1) = beta;
    
    figure('Name','Density Matrix Per Angle');
    pos = get(gcf,'pos');
    set(gcf,'pos',[pos(1) pos(2) 1400 1000]);
    figHand = gcf;
    % If no numOfSteps input
    if ~exist('numOfSteps','var')
        dens = plotDensityMatrix(polarization,endAngle);
    else
        % verify legal number of steps
        assert(numOfSteps > 0,'Number of steps must be positive number');

        % Angle vector definition
        dTheta = endAngle / numOfSteps;
        angles = 0:dTheta:endAngle;
        i = 1;
        key = ' ';
        while key~= "q"
            a = annotation('textbox',[0.1 0.1 0.14 0.04],'string',['Current angle: ' num2str(rad2deg(angles(i)))]);
            a.FontSize = 20;
            dens = plotDensityMatrix(polarization,angles(i));
            pause;
            key = get(figHand,'CurrentCharacter');
            if key == "o"
                if i > 1
                    i = i - 1;
                end
            else
                if key == "p"
                    if i < numel(angles)
                        i = i + 1;
                    end
                end
            end
            delete(a);
        end
        close(figHand);
    end
end

%% Plot Density Matrix
function mat = plotDensityMatrix(polarVec,theta)   
    % Define the matrices of the HWPs per this angle, the output
    % polarization vector, and the density matrix.
    % Calcite Operator Matrix:
    U_cal = createCalMat();
    U_cal_perp = createPerpCalMat();
    % HWP Operator Matrix:
    U_hwp_last = createHWPMatrix(-theta);
    U_hwp = createHWPMatrix(theta);
    % Output Polarization and density Matrix:
    outputPolarization = U_hwp_last*U_cal_perp*U_hwp*U_cal*polarVec;
    density = outputPolarization*(outputPolarization)';

    %% Plotting
    % Plot real part
    subplot(2,2,1);
    bar3(real(density));
    zlim([-1 1]);
    title('Real','FontSize',25);
    xticklabels({'<h0|','<h\tau|','<h2\tau|','<v0|','<v\tau|','<v2\tau|'});
    yticklabels({'|h0>','|h\tau>','|h2\tau>','|v0>','|v\tau>','|v2\tau>'});
    % Plot imaginary part
    subplot(2,2,2);
    bar3(imag(density));
    zlim([-1 1]);
    title('Imaginary','FontSize',25);
    xticklabels({'<h0|','<h\tau|','<h2\tau|','<v0|','<v\tau|','<v2\tau|'});
    yticklabels({'|h0>','|h\tau>','|h2\tau>','|v0>','|v\tau>','|v2\tau>'});
    % Plot absolute values
    subplot(2,2,[3,4]);
    bar3(abs(density));
    zlim([-1 1]);
    title('Absolute','FontSize',25);
    xticklabels({'<h0|','<h\tau|','<h2\tau|','<v0|','<v\tau|','<v2\tau|'});
    yticklabels({'|h0>','|h\tau>','|h2\tau>','|v0>','|v\tau>','|v2\tau>'});
    refresh;
    % Return Values
    mat = density;
end

% Create a HWP Operator Matrix in the angle theta
function U_hwp_mat = createHWPMatrix(theta)
    U_hwp_mat = zeros(2);
    U_hwp_mat(1,1) = cos(2*theta);
    U_hwp_mat(1,2) = sin(2*theta);
    U_hwp_mat(2,1) = sin(2*theta);
    U_hwp_mat(2,2) = -cos(2*theta);
    U_hwp_mat = kron(U_hwp_mat,eye(3));
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