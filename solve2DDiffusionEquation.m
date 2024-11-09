%--------------------------------------------------------------------------
%
%    solve2DDiffusionEquation.m
%
%    - A script to set up and solve the 2D diffusion equation for conduction
%    in a plate, with fixed temperature boundary conditions applied on the 
%    sides of the plate
%
%    [Chapter 2]
%
%    Author: Dr. Aidan Wimshurst (FluidMechanics101)
%
%    - Email: FluidMechanics101@gmail.com
%    - Web  : https://www.fluidmechanics101.com
%    - YouTube: Fluid Mechanics 101
%
%    Version 1.0.0 21/08/2021
%
%--------------------------------------------------------------------------
close all; clear; clc
%--------------------------------------------------------------------------
% User Inputs
%--------------------------------------------------------------------------
% Thermal Conductivity of the plate (W/mK)
conductivity = 100;

% Thickness of the plate (m)
thickness = 0.1;

% Length of the plate (x)
plateLength = 4.0;

% Width of the plate (y)
plateWidth = 4.0;

% Number of cells along the length (Nx)
numCellsLength = 4;

% Number of cells along the width (Ny)
numCellsWidth = 4;

% Temperature at the left end of the plate
temperatureLeft = 100;

% Temperature at the bottom of the plate
temperatureBottom = 150;

% Temperature at the right of the plate
temperatureRight = 200;

% Temperature at the top of the plate
temperatureTop = 250;

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
fprintf('   solve2DDiffusionEquation.m\n')
fprintf('\n')
fprintf(' - Fluid Mechanics 101\n')
fprintf(' - Author: Dr. Aidan Wimshurst\n')
fprintf(' - Contact: FluidMechanics101@gmail.com\n')
fprintf('\n')
fprintf(' [Chapter 2]\n')
fprintf('================================================\n')

%--------------------------------------------------------------------------
% 2. Create the mesh of cells
%--------------------------------------------------------------------------
fprintf(' Creating Mesh ...\n');
fprintf('------------------------------------------------\n');

% Start by calculating the number of cells
numCells = numCellsLength*numCellsWidth;

% Calculate the number of faces
numFaces = (numCellsLength+1)*(numCellsWidth+1);

% Coodinates of the faces
% - Assemble as a single long vector, rather than a grid.
% - repmat and repelem are useful for repeating regular patterns. 
xFacesPattern = linspace(0, plateLength, numCellsLength+1);
yFacesPattern = linspace(0, plateWidth, numCellsWidth+1);
xFaces = repmat(xFacesPattern, 1, numCellsWidth+1);
yFaces = repelem(yFacesPattern, numCellsLength+1);

% Coordinates of the centroids
xCentroidPattern = 0.5*(xFacesPattern(2:end) + xFacesPattern(1:end-1));
yCentroidPattern = 0.5*(yFacesPattern(2:end) + yFacesPattern(1:end-1));
xCentroids = repmat(xCentroidPattern , 1, numCellsWidth);
yCentroids = repelem(yCentroidPattern, numCellsLength);

% Distance between the cell centroids and the boundaries
dLeftBoundary = 2.0*(xCentroids(1) - xFacesPattern(1));
dTopBoundary = 2.0*(yCentroids(1) - yFacesPattern(1));
dRightBoundary = 2.0*(xFacesPattern(end) - xCentroids(end));
dBottomBoundary = 2.0*(yFacesPattern(end) - yCentroids(end));

% Assemble the distance vectors
dLeftPattern = [dLeftBoundary, xFacesPattern(2:end) - ...
    xFacesPattern(1:end-1)];
dRightPattern = [xFacesPattern(2:end) - xFacesPattern(1:end-1), ... 
     dRightBoundary];
dBottomPattern = [yFacesPattern(2:end)-yFacesPattern(1:end-1), ...
     dBottomBoundary];
dTopPattern = [dTopBoundary, yFacesPattern(2:end) - ...
     yFacesPattern(1:end-1)];
 
 dLeft = repmat(dLeftPattern(1:end-1), 1, numCellsWidth);
 dRight = repmat(dRightPattern(2:end), 1, numCellsWidth);
 dBottom = repelem(dBottomPattern(2:end), numCellsLength);
 dTop = repelem(dTopPattern(1:end-1), numCellsLength);
 
% Calculate the cell volumes. They are all have the same volume. 
cellLength = plateLength/numCellsLength;
cellWidth = plateWidth/numCellsWidth;
cellVolume = cellLength*cellWidth*thickness;

% Calculate the cross sectional area in the x and y directions
areaX = cellWidth*thickness;
areaY = cellLength*thickness;

% Identify the cells which have boundary faces. Give them an ID of 1. 
topBoundaryID = [ones(1, numCellsLength), ... 
    repmat(zeros(1, numCellsLength), 1, numCellsWidth - 1)];
bottomBoundaryID = [repmat(zeros(1, numCellsLength), 1, ...
    numCellsWidth - 1), ones(1, numCellsLength)];
leftBoundaryID = repmat([1.0, ...
    zeros(1, numCellsLength -1)], 1, numCellsWidth);
rightBoundaryID = repmat([zeros(1, numCellsLength -1), ...
    1.0], 1, numCellsWidth);

%--------------------------------------------------------------------------
% 3. Calculate matrix coefficients
%--------------------------------------------------------------------------
fprintf(' Calculating matrix coefficients ...\n');
fprintf('------------------------------------------------\n');

% Diffusive flux per unit area
DA_left = conductivity*areaX ./dLeft;
DA_right = conductivity*areaX ./dRight;
DA_bottom = conductivity*areaY ./dBottom;
DA_top = conductivity*areaY./dTop;

% Source term Su
% --------------
% The volume heat flux is the same for interior and boundary cells
Su = cellVolume*ones(1, numCells)*heatSourcePerVol;

% Add the contribution from each of the boundary faces
Su = Su + (2.0*temperatureLeft*leftBoundaryID.*DA_left);
Su = Su + (2.0*temperatureRight*rightBoundaryID.*DA_right);
Su = Su + (2.0*temperatureBottom*bottomBoundaryID.*DA_bottom);
Su = Su + (2.0*temperatureTop*topBoundaryID.*DA_top);

% Source term Sp
% ---------------
% The source term is zero for interior cells
Sp = zeros(1, numCells);

% Add the contribution from each of the boundary faces
Sp = Sp + (-2.0*DA_left.*leftBoundaryID);
Sp = Sp + (-2.0*DA_right.*rightBoundaryID);
Sp = Sp + (-2.0*DA_bottom.*bottomBoundaryID);
Sp = Sp + (-2.0*DA_top.*topBoundaryID);

% aL, aR, aT, aB
% --------------

% Only add contributions for interior cells
aL = DA_left.*(1- leftBoundaryID);
aR = DA_right.*(1- rightBoundaryID);
aB = DA_bottom.*(1- bottomBoundaryID);
aT = DA_top.*(1- topBoundaryID);

% Calculate ap from the other coefficients
aP = aL + aR + aB + aT - Sp;

%--------------------------------------------------------------------------
% 4. Print the setup
%--------------------------------------------------------------------------
if strcmp(printSetup, 'true')
    fprintf(' Summary: Set Up\n')
    fprintf('------------------------------------------------\n')
    fprintf('Cell | aL | aR  | aB  | aT  | Sp  |  Su  |  aP \n')
    fprintf('------------------------------------------------\n')
    for i = 1:numCells
       fprintf('%4i %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %6.1f\n', i, ...
           aL(i), aR(i), aB(i), aT(i), Sp(i), Su(i), aP(i)); 
    end
    fprintf('------------------------------------------------\n');
end

%--------------------------------------------------------------------------
% 5.0 Create the matrices
%--------------------------------------------------------------------------
fprintf(' Assembling Matrices\n');
fprintf('------------------------------------------------\n');

% Start by creating an empty A matrix and an empty B matrix
Amatrix = zeros(numCells, numCells);
BVector = zeros(numCells, 1);

% Populate the matrix one row at a time (i.e one cell at a time)
%
% NOTE: this method is deliberately inefficient for this problem
% but it is useful for learning purposes. We could populate the
% diagonals and the off-diagonals directly.

for i = 1:numCells

    % Do the A matrix diagonal coefficients
    Amatrix(i, i) = aP(i);

    % Do the B matrix coefficients
    BVector(i) = Su(i);

    % Does the cell have a left boundary?
    if (leftBoundaryID(i) == 0.0)
        Amatrix(i, i-1) = -1.0*aL(i);
    end

    % Does the cell have a right boundary?
    if (rightBoundaryID(i) == 0.0)
        Amatrix(i, i+1) = -1.0*aR(i);
    end

    % Does the cell have a bottom boundary?
    if (bottomBoundaryID(i) == 0.0)
        Amatrix(i, i+numCellsLength) = -1.0*aB(i);
    end

    % Does the cell have a top boundary?
    if (topBoundaryID(i) == 0.0)
        Amatrix(i, i-numCellsLength) = -1.0*aT(i);
    end
end

%--------------------------------------------------------------------------
% 5.0 Print the matrices
%--------------------------------------------------------------------------

if strcmp(printSetup, 'true')
    fprintf('A matrix:\n')
    for i = 1:numCells
        for j = 1:numCells
           fprintf('%6.1f ', Amatrix(i, j)); 
        end
        fprintf('\n')
    end
    fprintf('B vector:\n')
    for i = 1:numCells
       fprintf('%6.1f\n', BVector(i)); 
    end
    fprintf('\n')
    fprintf('------------------------------------------------\n')
end

%--------------------------------------------------------------------------
% 6.0 Solve the matrices
%--------------------------------------------------------------------------

fprintf(' Solving ...\n')
fprintf('------------------------------------------------\n')

% Use MATLAB's default linear algebra solver (AX = B)
Tvector = Amatrix \ BVector;

fprintf(' Equations Solved.\n')
fprintf('------------------------------------------------\n')

%--------------------------------------------------------------------------
% 7.0 Reshape the temperature vector into a grid
%--------------------------------------------------------------------------
Tgrid = reshape(Tvector, [numCellsLength, numCellsWidth])';

%--------------------------------------------------------------------------
% 8.0 Print the Solution
%--------------------------------------------------------------------------
if strcmp(printSolution, 'true')
    fprintf(' Solution: Temperature Vector\n');
    fprintf('------------------------------------------------\n');
    fprintf('T vector:\n')
    for i = 1:numCells
       fprintf('%10.6f\n', Tvector(i)); 
    end
    fprintf('\n')
    fprintf('------------------------------------------------\n')
end
%--------------------------------------------------------------------------
% 8.0 Heat Fluxes
%--------------------------------------------------------------------------

% Calculate the temperature differences
% - To do this we need to stack on the boundary temperatures onto the grid
Tleftrightshift = [temperatureLeft*ones(numCellsWidth, 1), ...
    Tgrid, temperatureRight*ones(numCellsWidth,1)]; 
Ttopbottomshift = vertcat(temperatureTop*ones(1, numCellsLength), ...
    Tgrid, temperatureBottom*ones(1, numCellsLength));

% Now we can calculate the temperature differences
deltaTleft = Tleftrightshift(:,2:end-1)-Tleftrightshift(:,1:end-2);
deltaTright = Tleftrightshift(:,3:end)-Tleftrightshift(:,2:end-1);
deltaTtop = Ttopbottomshift(1:end-2,:)-Ttopbottomshift(2:end-1,:);
deltaTbottom = Ttopbottomshift(2:end-1,:) - Ttopbottomshift(3:end,:);

% We now need to calculate the diffusive heat flux (DA) on each face
% - Start by reshaping the DA vectors into a grid of the correct size
DA_left_grid = reshape(DA_left, [numCellsWidth, numCellsLength]);
DA_right_grid = reshape(DA_right, [numCellsWidth, numCellsLength]);
DA_top_grid = reshape(DA_top, [numCellsWidth, numCellsLength]);
DA_bottom_grid = reshape(DA_bottom, [numCellsWidth, numCellsLength]);

% Calculate the boundary face fluxes
DA_left_boundary = (2.0*conductivity*areaX/dLeftBoundary)* ...
    ones(numCellsWidth,1);
DA_right_boundary = (2.0*conductivity*areaX/dRightBoundary)*...
    ones(numCellsWidth,1);
DA_top_boundary = (2.0*conductivity*areaY/dTopBoundary)*...
    ones(1, numCellsLength);
DA_bottom_boundary = (2.0*conductivity*areaY/dBottomBoundary)*...
    ones(1, numCellsLength);

% Now stack on the boundary face fluxes to the grid
DA_left_shift = [DA_left_boundary, DA_left_grid(:,2:end)];
DA_right_shift = [DA_right_grid(:,1:end-1), DA_right_boundary];
DA_top_shift = vertcat(DA_top_boundary, DA_top_grid(2:end,:));
DA_bottom_shift = vertcat(DA_bottom_grid(1:end-1,:), DA_top_boundary);

% Unit normal vectors
normalsLeftGrid = -1.0*ones(numCellsWidth, numCellsLength);
normalsRightGrid = ones(numCellsWidth, numCellsLength);
normalsBottomGrid = -1.0*ones(numCellsWidth, numCellsLength);
normalsTopGrid = ones(numCellsWidth, numCellsLength);

% Now we can compute the heat fluxes
heatFluxLeft = -DA_left_shift.* deltaTleft.*normalsLeftGrid;
heatFluxRight = -DA_right_shift.* deltaTright.*normalsRightGrid;
heatFluxTop = -DA_top_shift.* deltaTtop.*normalsTopGrid;
heatFluxBottom = -DA_bottom_shift.* deltaTbottom.*normalsBottomGrid;

% Calculate the volumetric heat source in each cell
sourceVol = heatSourcePerVol*cellVolume*ones(numCellsWidth, ...
    numCellsLength);

% Calculate the error in the heat flux balance in each cell
error = (sourceVol - heatFluxLeft - heatFluxRight - heatFluxTop - ...
    heatFluxBottom);

% Reshape the matrices into vectors for printing
heatFluxLeftVector = reshape(heatFluxLeft', 1, []);
heatFluxRightVector = reshape(heatFluxRight', 1, []);
heatFluxTopVector = reshape(heatFluxTop', 1, []);
heatFluxBottomVector = reshape(heatFluxBottom', 1, []);
sourceVolVector = reshape(sourceVol', 1, []);
errorVector = reshape(error', 1, []);

%--------------------------------------------------------------------------
% 9.0 Print the Heat Flux Balance
%--------------------------------------------------------------------------

if strcmp(printSolution, 'true')

    fprintf(' Heat Fluxes\n');
    fprintf('---------------------------------------------------\n');
    fprintf('Cell |  QL   |  QR   |  QT   |  QB   |  SV   |  Err\n');
    fprintf('---------------------------------------------------\n');
    for i = 1:numCells
        fprintf('%4i %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f\n',  i, ...
            heatFluxLeftVector(i), heatFluxRightVector(i), ...
            heatFluxTopVector(i), heatFluxBottomVector(i), ...
            sourceVolVector(i), errorVector(i));
    end
    fprintf('---------------------------------------------------\n');
end

% Sum the heat fluxes across the boundary faces to give the total heat flux
% across each boundary
heatFluxLeftBoundaryTotal = sum(leftBoundaryID.*heatFluxLeftVector);
heatFluxRightBoundaryTotal = sum(rightBoundaryID.*heatFluxRightVector);
heatFluxBottomBoundaryTotal = sum(bottomBoundaryID.* heatFluxBottomVector);
heatFluxTopBoundaryTotal = sum(topBoundaryID.*heatFluxTopVector);
heatFluxBoundaryTotal = (heatFluxLeftBoundaryTotal + ... 
                        heatFluxRightBoundaryTotal + ...
                        heatFluxTopBoundaryTotal + ...
                        heatFluxBottomBoundaryTotal);
heatGeneratedTotal = sum(sourceVolVector);

if strcmp(printSolution, 'true')

    fprintf(' Boundary Heat Flux Balance\n');
    fprintf('---------------------------------------------------\n');
    fprintf(' Left      : %7.1f [W]\n',  heatFluxLeftBoundaryTotal);
    fprintf(' Right     : %7.1f [W]\n',  heatFluxRightBoundaryTotal);
    fprintf(' Bottom    : %7.1f [W]\n',  heatFluxBottomBoundaryTotal);
    fprintf(' Top       : %7.1f [W]\n',  heatFluxTopBoundaryTotal);
    fprintf(' Total     : %7.1f [W]\n',  heatFluxBoundaryTotal);
    fprintf('---------------------------------------------------\n');
    fprintf(' Generated : %7.1f [W]\n',  heatGeneratedTotal);
    fprintf('---------------------------------------------------\n');
    fprintf(' Error     : %7.1f [W]\n',  (heatFluxBoundaryTotal - ...
        heatGeneratedTotal));
    fprintf('---------------------------------------------------\n');
end

%--------------------------------------------------------------------------
% 10.0 Interpolate onto node grid for plotting
%--------------------------------------------------------------------------

% Interpolate the solution on the interior nodes from the CFD solution
Tleftrightnodes = 0.5*(Tgrid(:,2:end)+Tgrid(:,1:end-1));
Tinternalnodes = 0.5*(Tleftrightnodes(2:end,:) + ...
    Tleftrightnodes(1:end-1,:));

% Interpolate the boundary temperatures in the corners
temperatureTopLeftCorner = 0.5*(temperatureTop+temperatureLeft);
temperatureTopRightCorner = 0.5*(temperatureTop+temperatureRight);
temperatureBottomLeftCorner = 0.5*(temperatureBottom+temperatureLeft);
temperatureBottomRightCorner = 0.5*(temperatureBottom+temperatureRight);

% Assemble the temperatures on all the boundary nodes
temperatureTopVector = [temperatureTopLeftCorner, ...
    temperatureTop*ones(1, numCellsLength-1), temperatureTopRightCorner];
temperatureBottomVector = [temperatureBottomLeftCorner, ... 
    temperatureBottom*ones(1, numCellsLength-1), ...
    temperatureBottomRightCorner];
temperatureLeftVector = temperatureLeft*ones(numCellsWidth-1, 1);
temperatureRightVector = temperatureRight*ones(numCellsWidth-1,1);

% Assemble the temperature on all of the nodes together as one grid
Tnodes = vertcat(temperatureTopVector, [temperatureLeftVector, ...
    Tinternalnodes, temperatureRightVector], temperatureBottomVector);

% X and Y coordinates of the nodes
xNodes = reshape(xFaces, [numCellsLength+1, numCellsWidth+1])';
yNodes = flipud(reshape(yFaces, [numCellsLength+1, numCellsWidth+1])');

%--------------------------------------------------------------------------
% 11.0 Plotting
%--------------------------------------------------------------------------

% Plot the data if desired
if strcmp(plotOutput,'true')

    fprintf(' Plotting ...\n')
    fprintf('------------------------------------------------\n')
    
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
    
    fig1 = figure('Name', '2D Diffusion');
    box on;
    colormap(flipud(autumn));
    contourf(xNodes, yNodes, Tnodes, 'k-o', ...
        'linewidth', lineWidth + 0.2);
    xlabel('x [m]', 'FontSize', fontSize)
    ylabel('y [m]', 'FontSize', fontSize)
    c = colorbar('eastoutside');
    c.Label.String = 'T [^{\circ} C]';
    caxis([100, 260])
    c.Ticks = [100, 120, 140, 160, 180, 200, 220, 240, 260];
    set(gca, 'fontsize', fontSize);
    set(gca, 'linewidth', lineWidth);
    set(gcf, 'position', [x0, y0, figSizeXpixels, figSizeYpixels]);
    set(gca,'XTick',0:1:4)
    
end

%==========================================================================
% END OF FILE
%==========================================================================


