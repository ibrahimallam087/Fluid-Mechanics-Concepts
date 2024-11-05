###############################################################################
#
#    neumannBoundaryConditions.py
#
#    - A script to set up and solve the 1D diffusion equation for conduction
#    in a bar, with heat flux boundary conditions at the left end and
#    fixed temperature conditions at the right end
#
#    [Chapter 1]
#
#    Author: Dr. Aidan Wimshurst (FluidMechanics101)
#
#    - Email: FluidMechanics101@gmail.com
#    - Web  : https://www.fluidmechanics101.com
#    - YouTube: Fluid Mechanics 101
#
#    Version 2.0.0 15/07/2019
#
#    Version 2.0.3 15/04/2022
#    - Updated for python 2 & python 3 compatability.
#
###############################################################################
#
# Required Modules for Matrices and Plotting
#
###############################################################################
import numpy as np
import matplotlib.pyplot as plt
###############################################################################
#
#    User Inputs
#
###############################################################################

# Thermal Conductivity of the bar (W/mK)
cond = 100

# Cross-sectional Area of the bar (m2)
area = 0.1

# Length of the bar (m)
barLength = 5

# Number of cells in the mesh
nCells = 5

# Heat flux from the left end of the bar (W/m2)
heatFluxLeftBoundary = 100

# Temperature at the right hand side of the bar (deg C)
tempRight = 200

# Heat source per unit volume (W/m3)
heatSourcePerVol = 1000

# Plot the data?
plotOutput = 'true'

# Print the set up data? (table of coefficients and matrices)
printSetup = 'true'

# Print the solution output (the temperature vector)
printSolution = 'true'

# =============================================================================
# Code Begins Here
# =============================================================================

# Define a line for printing (to ensure they have the same length)
lineSingle = '------------------------------------------------'
lineDouble = '================================================'

print (lineDouble)
print ('')
print ('   neumannBoundaryConditions.py')
print ('')
print (' - Fluid Mechanics 101')
print (' - Author: Dr. Aidan Wimshurst')
print (' - Contact: FluidMechanics101@gmail.com')
print ('')
print (' [Chapter 1]')
print (lineDouble)

###############################################################################
#
#     Create the mesh of cells
#
###############################################################################

print (' Creating Mesh')
print (lineSingle)

# Start by calculating the coordinates of the cell faces
xFaces = np.linspace(0, barLength, nCells+1)

# Calculate the coordinates of the cell centroids
xCentroids = 0.5*(xFaces[1:] + xFaces[0:-1])

# Calculate the length of each cell
cellLength = xFaces[1:] - xFaces[0:-1]

# Calculate the distance between cell centroids
dCentroids = xCentroids[1:] - xCentroids[0:-1]

# For the boundary cell on the left, the distance dLP is double the distance
# from the cell centroid to the boundary face
dLeft = 2*(xCentroids[0] - xFaces[0])

# For the boundary cell on the right, the distance dPR is double the distance
# from the cell centroid to the boundary cell face
dRight = 2*(xFaces[-1] - xCentroids[-1])

# Append these to the vector of distances
dCentroids = np.hstack([dLeft, dCentroids, dRight])

# Distance between the centroid and the centroid on the left
dCentroidLeft = dCentroids[0:-1]

# Distance between the centroid and the centroid on the right
dCentroidRight = dCentroids[1:]

# Calculate the area of the left faces
areaLeftFaces = area*np.ones(nCells)

# Calculate the area of the right faces
areaRightFaces = area*np.ones(nCells)

# Compute the cell volume
cellVolumes = cellLength*area

###############################################################################
#
#     Material Properties
#
###############################################################################

# Assign the thermal conductivity to all of the cell faces
conductivityFaces = cond*np.ones(nCells+1)

# Extract a vector for the conductivity of the left faces
conductivityLeftFaces = conductivityFaces[0:-1]

# Extract a vector for the conductivity of the right faces
conductivityRightFaces = conductivityFaces[1:]

###############################################################################
#
#    Calculate the Matrix Coefficients
#
###############################################################################

print (' Calculating Matrix Coefficients')
print (lineSingle)

# Diffusive flux per unit area
DA_LeftFaces = np.divide(
    np.multiply(conductivityLeftFaces, areaLeftFaces), dCentroidLeft)
DA_RightFaces = np.divide(
    np.multiply(conductivityRightFaces, areaRightFaces), dCentroidRight)

# Calculate the source term Su
Su = heatSourcePerVol*cellVolumes

# Assign sources to the left and right hand boundaries
Su[0] = Su[0] - heatFluxLeftBoundary*area
Su[-1] = Su[-1] + 2.0*np.copy(DA_RightFaces[-1])*tempRight

# Calculate the source term Sp
Sp = np.zeros(nCells)

# Assign sources to the right hand boundary
Sp[0] = 0.0
Sp[-1] = -2.0*np.copy(DA_RightFaces[-1])

# Left Coefficient (aL)
aL = np.copy(DA_LeftFaces)

# Right Coefficient (aR)
aR = np.copy(DA_RightFaces)

# Set the first element of aL to zero. It is a boundary face
aL[0] = 0

# Set the last element of aR to zero. It is a boundary face
aR[-1] = 0

# Create the central coefficients
aP = np.copy(aL) + np.copy(aR) - np.copy(Sp)

###############################################################################
#
#     Create the matrices
#
###############################################################################

print (' Assembling Matrices')
print (lineSingle)

# Start by creating an empty A matrix and an empty B matrix
Amatrix = np.zeros([nCells, nCells])
BVector = np.zeros(nCells)

# Populate the matrix one row at a time (i.e one cell at a time)

#
# NOTE: this method is deliberately inefficient for this problem
#        but it is useful for learning purposes. We could populate the
#        diagonals and the off-diagonals directly.

for i in range(nCells):

    # Do the A matrix Coefficients

    # Left boundary Cell
    if (i == 0):

        Amatrix[i, i] = aP[i]
        Amatrix[i, i+1] = -1.0*aR[i]

    # Right Boundary Cell
    elif(i == nCells-1):

        Amatrix[i, i-1] = -1.0*aL[i]
        Amatrix[i, i] = aP[i]

    # Interior Cells
    else:

        Amatrix[i, i-1] = -1.0*aL[i]
        Amatrix[i, i] = aP[i]
        Amatrix[i, i+1] = -1.0*aR[i]

    # Do the B matrix coefficients
    BVector[i] = Su[i]

###############################################################################
#
#     Print the setup
#
###############################################################################

if (printSetup == 'true'):

    print (' Summary: Set Up')
    print (lineSingle)
    print ('Cell | aL | aR  | ap  |  Sp  |   Su ')
    print (lineSingle)
    for i in range(nCells):
        print (('%4i %5.1f %5.1f %5.1f %5.1f %7.1f ' % (
            i+1, aL[i], aR[i], aP[i], Sp[i], Su[i])))
    print (lineSingle)
    np.set_printoptions(linewidth=np.inf)
    print ('A matrix:')
    print (lineSingle)
    print (Amatrix)
    print ('B vector')
    print (lineSingle)
    print (BVector)
    print (lineSingle)

###############################################################################
#
#     Solve the matrices
#
###############################################################################

print (' Solving ...')
print (lineSingle)

# Use the built-in python solution module
Tvector = np.linalg.solve(Amatrix, BVector)

print (' Equations Solved')
print (lineSingle)

###############################################################################
#
#     Print the Results
#
###############################################################################

if (printSolution == 'true'):

    print (' Solution: Temperature Vector')
    print (lineSingle)
    print (Tvector)
    print (lineSingle)

###############################################################################
#
#     Heat Fluxes
#
###############################################################################

# Use central differencing to calculate the left boundary temperature
tempLeft = (Tvector[0]
            - (heatFluxLeftBoundary*area)/(2.0*np.copy(DA_LeftFaces[0])))

# Stack on the boundary temperatures to the temperature vector
temperatureStack = np.hstack([tempLeft, np.copy(Tvector), tempRight])

# Calculate the temperature differences
temperatureDifferencesLeft = temperatureStack[1:-1] - temperatureStack[0:-2]
temperatureDifferencesRight = temperatureStack[2:] - temperatureStack[1:-1]

# Unit normal vectors
normalsLeft = -1.0*np.ones(nCells)
normalsRight = np.ones(nCells)

# Calculate the heat fluxes
heatFluxLeft = -1.0*np.prod(
        [normalsLeft, temperatureDifferencesLeft, DA_LeftFaces], 0)
heatFluxRight = -1.0*np.prod(
        [normalsRight, temperatureDifferencesRight, DA_RightFaces], 0)

# For the left and right boundary face the heat flux
# is 2*DA*(temperatureDifference), as the distance is halved.
heatFluxLeft[0] *= 2.0
heatFluxRight[-1] *= 2.0

# Calculate the heat source in each cell
heatSource = heatSourcePerVol*cellVolumes*np.ones(nCells)

# Calculate the heat flux error in each cell
heatBalanceError = heatSource - heatFluxLeft - heatFluxRight

if (printSolution == 'true'):

    print (' Heat Fluxes')
    print (lineSingle)
    print ('Cell |  QL   |  QR   |  SV   |  Err')
    print (lineSingle)
    for i in range(nCells):
        print (('%4i %7.1f %7.1f %7.1f %7.1f' % (
            i+1, heatFluxLeft[i], heatFluxRight[i],
            heatSource[i], heatBalanceError[i])))
    print (lineSingle)

###############################################################################
#
#     Plot the solution
#
###############################################################################

# Plot the data if desired
if (plotOutput == 'true'):

    print (' Plotting ...')
    print (lineSingle)

    # Append the boundary temperature values to the vector, so we can
    # plot the complete solution
    xPlotting = np.hstack([xFaces[0], np.copy(xCentroids), xFaces[-1]])

    # Assemble the analytical solution for comparison
    xAnalytical = np.linspace(0, barLength, 100)
    temperatureAnalytical = ((heatSourcePerVol/(2.0*cond))
                             * (barLength*barLength*np.ones(len(xAnalytical))
                             - np.square(xAnalytical))
                             + (heatFluxLeftBoundary/cond)*(xAnalytical
                             - barLength*np.ones(len(xAnalytical)))
                             + tempRight)

    # Configure the plot to look how you want
    # - The plot is currently 3.1 inches by 2.5 inches, so that 2 images can
    #     fit side by side on a page.
    fontSize = 11
    fontSizeLegend = 11
    lineWidth = 1.5
    tickPad = 8
    tickPad2 = 16
    labelPadY = 3
    labelPadX = 2
    boxPad = 2.5
    tickLength = 4
    markerSize = 4

    # Colours - Can use rgb or html
    lightBlue = '#bfc8d1'
    shadeBlue = '#8091a4'
    darkBlue = '#002147'

    # Use latex 'CMU sans-serif' font in the plots.
    plt.rc('font', family='serif')
    plt.rcParams['axes.linewidth'] = lineWidth
    plt.rcParams["figure.figsize"] = (3.1, 2.5)

    fig1 = plt.figure(1)
    ax = fig1.add_subplot(111)
    fig1.tight_layout(pad=boxPad)
    ax.plot(xPlotting, temperatureStack, '-o',
            linewidth=lineWidth, label='CFD', color=darkBlue,
            markersize=markerSize)
    ax.plot(xAnalytical, temperatureAnalytical, '--',
            linewidth=lineWidth, label='Analytical', color=shadeBlue)
    plt.xlabel(r'$x$ [m]', fontsize=fontSize, labelpad=labelPadX)
    plt.ylabel(r'$T$ [$^{\circ}$ C]', fontsize=fontSize, labelpad=labelPadY)
    plt.yticks(fontsize=fontSize)
    plt.xticks(fontsize=fontSize)
    plt.xlim([xFaces[0], xFaces[-1]])
    plt.ylim([200, 350])
    leg = plt.legend(fontsize=fontSizeLegend, loc='best', fancybox=False)
    leg.get_frame().set_linewidth(lineWidth)
    ax.tick_params(which='both', direction='in', length=6,
                   width=lineWidth, gridOn=False, pad=tickPad, color=darkBlue)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.spines['bottom'].set_color(darkBlue)
    ax.spines['top'].set_color(darkBlue)
    ax.spines['right'].set_color(darkBlue)
    ax.spines['left'].set_color(darkBlue)
    plt.show()

###############################################################################
#
#     End of File
#
###############################################################################
