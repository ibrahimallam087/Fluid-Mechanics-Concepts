%--------------------------------------------------------------------------
%
%    neumannBoundaryConditions.m
%
%    - A script to set up and solve the 1D diffusion equation for 
%    conduction in a bar, with heat flux boundary conditions at the left
%    end and fixed temperature conditions at the right end.
%
%    [Chapter 1]
%
%    Author: Dr. Aidan Wimshurst (FluidMechanics101)
%
%    - Email: FluidMechanics101@gmail.com
%    - Web  : https://www.fluidmechanics101.com
%    - YouTube: Fluid Mechanics 101
%
%    Version 1.0.0 18/08/2021
%
%--------------------------------------------------------------------------
close all; clear; clc
%--------------------------------------------------------------------------
% User Inputs
%--------------------------------------------------------------------------

% Thermal Conductivity of the bar (W/mK)
cond = 100;

% Cross-sectional Area of the bar (m2)
area = 0.1;

% Length of the bar (m)
barLength = 5;

% Number of cells in the mesh
nCells = 5;

% Heat flux from the left end of the bar (W/m2)
heatFluxLeftBoundary = 100;

% Temperature at the right hand side of the bar (deg C)
tempRight = 200;

% Heat source per unit volume (W/m3)
heatSourcePerVol = 1000;

% Plot the data?
plotOutput = 'true';

% Print the set up data? (table of coefficients and matrices)
printSetup = 'true';

% Print the solution output (the temperature vector)
printSolution = 'true';

%==========================================================================
% Code Begins Here
%==========================================================================
%--------------------------------------------------------------------------
% 1. Print messages
%--------------------------------------------------------------------------

fprintf('================================================\n');
fprintf('\n');
fprintf('   solve1DDiffusionEquation.m\n')
fprintf('\n')
fprintf(' - Fluid Mechanics 101\n')
fprintf(' - Author: Dr. Aidan Wimshurst\n')
fprintf(' - Contact: FluidMechanics101@gmail.com\n')
fprintf('\n')
fprintf(' [Exercise 1: Chapter 2]\n')
fprintf('================================================\n')

%--------------------------------------------------------------------------
% 2. Create the mesh of cells
%--------------------------------------------------------------------------
fprintf(' Creating Mesh ...\n');
fprintf('------------------------------------------------\n');

% Start by calculating the coordinates of the cell faces
xFaces = linspace(0, barLength, nCells+1);

% Calculate the coordinates of the cell centroids
xCentroids = 0.5*(xFaces(2:end) + xFaces(1:end-1));

% Calculate the length of each cell
cellLength = xFaces(2:end) - xFaces(1:end-1);

% Calculate the distance between cell centroids
dCentroids = xCentroids(2:end) - xCentroids(1:end-1);

% For the boundary cell on the left, the distance is double the distance
% from the cell centroid to the boundary face
dLeft = 2*(xCentroids(1) - xFaces(1));

% For the boundary cell on the right, the distance is double the distance from
% the cell centroid to the boundary cell face
dRight = 2*(xFaces(end) - xCentroids(end));

% Append these to the vector of distances
dCentroids = [dLeft, dCentroids, dRight];

% Distance between the centroid and the centroid on the left
dCentroidLeft = dCentroids(1:end-1);

% Distance between the centroid and the centroid on the right
dCentroidRight = dCentroids(2:end);

% Calculate the area of the left faces
areaLeftFaces = area*ones(1, nCells);

% Calculate the area of the right faces
areaRightFaces = area*ones(1, nCells);

% Compute the cell volume
cellVolumes = cellLength*area;

%--------------------------------------------------------------------------
% 3.0 Material Properties
%--------------------------------------------------------------------------
% Assign the thermal conductivity to all of the cell faces
conductivityFaces = cond*ones(nCells+1);

% Extract a vector for the conductivity of the left faces
conductivityLeftFaces = conductivityFaces(1, 1:end-1);

% Extract a vector for the conductivity of the right faces
conductivityRightFaces = conductivityFaces(1, 2:end);

%--------------------------------------------------------------------------
% 4.0 Calculate the Matrix Coefficients
%--------------------------------------------------------------------------
fprintf(' Calculating Matrix Coefficients\n');
fprintf('------------------------------------------------\n');

% Diffusive flux per unit area
DA_LeftFaces = conductivityLeftFaces .* areaLeftFaces ./ dCentroidLeft;
DA_RightFaces = conductivityRightFaces .* areaRightFaces ./ dCentroidRight;

% Calculate the source term Su
Su = heatSourcePerVol*cellVolumes;

% Assign sources to the left and right hand boundaries
Su(1) = Su(1) - heatFluxLeftBoundary*area;
Su(end) = Su(end) + 2.0*DA_RightFaces(end)*tempRight;

% Calculate the source term Sp
Sp = zeros(1, nCells);

% Assign sources to the right hand boundary
Sp(1) = 0.0;
Sp(end) = -2.0*DA_RightFaces(end);

% Left Coefficient (aL)
aL = DA_LeftFaces;

% Right Coefficient (aR)
aR = DA_RightFaces;

% Set the first element of aL to zero. It is a boundary face
aL(1) = 0;

% Set the last element of aR to zero. It is a boundary face
aR(end) = 0;

% Create the central coefficients
aP = aL + aR - Sp;

%--------------------------------------------------------------------------
% 5.0 Create the matrices
%--------------------------------------------------------------------------

fprintf(' Assembling Matrices\n');
fprintf('------------------------------------------------\n');

% Start by creating an empty A matrix and an empty B matrix
Amatrix = zeros(nCells, nCells);
BVector = zeros(nCells, 1);

% Populate the matrix one row at a time (i.e one cell at a time)
%
% NOTE: this method is deliberately inefficient for this problem
% but it is useful for learning purposes. We could populate the
% diagonals and the off-diagonals directly.

for i = 1:nCells

    % Do the A matrix Coefficients

    % Left boundary Cell
    if (i == 1)

        Amatrix(i, i) = aP(i);
        Amatrix(i, i+1) = -1.0*aR(i);

    % Right Boundary Cell
    elseif(i == nCells)

        Amatrix(i, i-1) = -1.0*aL(i);
        Amatrix(i, i) = aP(i);

    % Interior Cells
    else

        Amatrix(i, i-1) = -1.0*aL(i);
        Amatrix(i, i) = aP(i);
        Amatrix(i, i+1) = -1.0*aR(i);
    end
    
    % Do the B matrix coefficients
    BVector(i) = Su(i);
end

%--------------------------------------------------------------------------
% 4.0 Print the setup
%--------------------------------------------------------------------------
if strcmp(printSetup, 'true')
    fprintf(' Summary: Set Up\n')
    fprintf('------------------------------------------------\n')
    fprintf('Cell | aL | aR  | aP  |  Sp  |  Su \n')
    fprintf('------------------------------------------------\n')
    for i = 1:nCells
       fprintf('%4i %5.1f %5.1f %5.1f %5.1f %6.1f\n', i, aL(i), aR(i), aP(i), ...
           Sp(i), Su(i)); 
    end
    fprintf('------------------------------------------------\n')
end

%--------------------------------------------------------------------------
% 5.0 Solve the matrices
%--------------------------------------------------------------------------

fprintf(' Solving ...\n')
fprintf('------------------------------------------------\n')

% Use MATLAB's default linear algebra solver (AX = B)
Tvector = Amatrix \ BVector;

fprintf(' Equations Solved.\n')
fprintf('------------------------------------------------\n')

%--------------------------------------------------------------------------
% 6.0 Print the Solution
%--------------------------------------------------------------------------
if strcmp(printSolution, 'true')
    fprintf(' Solution: Temperature Vector\n');
    fprintf('------------------------------------------------\n');
    fprintf('T vector:\n')
    for i = 1:nCells
       fprintf('%10.6f\n', Tvector(i)); 
    end
    fprintf('\n')
    fprintf('------------------------------------------------\n')
end
%--------------------------------------------------------------------------
% 7.0 Heat Fluxes
%--------------------------------------------------------------------------
% Use central differencing to calculate the left boundary temperature
tempLeft = Tvector(1) - (heatFluxLeftBoundary.*area)./(2.0*DA_LeftFaces(1));

% Stack on the boundary temperatures to the temperature vector
temperatureStack = [tempLeft; Tvector; tempRight];

% Calculate the temperature differences
temperatureDifferencesLeft = temperatureStack(2:end-1) - temperatureStack(1:end-2);
temperatureDifferencesRight = temperatureStack(3:end) - temperatureStack(2:end-1);

% Unit normal vectors
normalsLeft = -1.0*ones(nCells, 1);
normalsRight = ones(nCells, 1);

% Calculate the heat fluxes
heatFluxLeft = -1.0.*normalsLeft.*temperatureDifferencesLeft.*DA_LeftFaces';
heatFluxRight = -1.0*normalsRight.*temperatureDifferencesRight.* ... 
    DA_RightFaces';

% For the left and right boundary face the heat flux 
% is 2*DA*(temperatureDifference), as the distance is halved. 
heatFluxLeft(1) = heatFluxLeft(1)*2.0;
heatFluxRight(end) = heatFluxRight(end)*2.0;

% Calculate the heat source in each cell
heatSource = heatSourcePerVol*cellVolumes;

% Calculate the heat flux error in each cell
heatBalanceError = heatSource - heatFluxLeft - heatFluxRight;

if strcmp(printSolution, 'true')

    fprintf(' Heat Fluxes\n');
    fprintf('------------------------------------------------\n');
    fprintf('Cell |  QL   |  QR   |  SV   |  Err\n');
    fprintf('------------------------------------------------\n');
    for i = 1:nCells
        fprintf('%4i %7.1f %7.1f %7.1f %7.1f\n',  i, ...
            heatFluxLeft(i), heatFluxRight(i), ... 
            heatSource(i), heatBalanceError(i));
    end
    fprintf('------------------------------------------------\n');
end

%--------------------------------------------------------------------------
% 8.0 Plot the solution
%--------------------------------------------------------------------------

% Plot the data if desired
if strcmp(plotOutput,'true')

    fprintf(' Plotting ...\n')
    fprintf('------------------------------------------------\n')
    
    % Append the boundary temperature values to the vector, so we can
    % plot the complete solution
    xPlotting = [xFaces(1), xCentroids, xFaces(end)];
    temperaturePlotting = [tempLeft; Tvector; tempRight];
    
    % Assemble the analytical solution for comparison
    xAnalytical = linspace(0, barLength, 100);                    
    temperatureAnalytical = ((heatSourcePerVol/(2.0*cond)) * ...
                            (barLength*barLength*ones(1, length(xAnalytical)) - ...
                            power(xAnalytical, 2)) + ... 
                            (heatFluxLeftBoundary/cond)*(xAnalytical - ... 
                            barLength*ones(1, length(xAnalytical))) + ...
                            tempRight);

    % Figure Size Parameters
    figSizeXcm = 8.5;
    aspectRatio = 4.0/3.0;
    figSizeYcm = figSizeXcm / aspectRatio;
    figSizeXpixels = figSizeXcm * 37.79527559055118;
    figSizeYpixels = figSizeYcm * 37.79527559055118;

    % Figure font, text and line widths
    fontSize = 10;
    fontChoice = 'Arial';
    lineWidth = 1.0;
    markerSize = 4;
    
    % Colours for the line plots
    colour1 = 'black';
    colour2 = [0.0, 0.129, 0.2784];

    % Adjust anchor for where figure shows up on screen (px)
    x0 = 100;
    y0 = 100;

    fig1 = figure('Name', '1D Diffusion');
    box on;
    hold on;
    plot(xPlotting, temperaturePlotting, 'k-o', 'Color', colour2, ...
        'MarkerSize', markerSize, 'linewidth', lineWidth + 0.2);
    plot(xAnalytical, temperatureAnalytical, 'k--', 'Color', colour2, ...
        'MarkerSize', markerSize, 'linewidth', lineWidth + 0.2);
    hold off;
    xlabel('x [m]', 'FontSize', fontSize)
    ylabel('T [^{\circ} C]', 'FontSize', fontSize)
    grid on;
    set(gca, 'fontsize', fontSize);
    set(gca, 'linewidth', lineWidth);
    set(gcf, 'position', [x0, y0, figSizeXpixels, figSizeYpixels]);
    legend({'CFD', 'Analytical'}, 'NumColumns', 1, 'Location', 'Best')
    
end

%==========================================================================
% END OF FILE
%==========================================================================
