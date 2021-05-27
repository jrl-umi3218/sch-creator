import numpy
from numpy.core.defchararray import add, index
from numpy.lib.function_base import angle, extract, flip
from numpy.lib.scimath import arccos
import yaml
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
from matplotlib.widgets import Slider
import numpy as np
from numpy import linalg as LA

########## functions ##########
def findAngle(vector1, vector2):
    uVector1 = vector1 / LA.norm(vector1)
    uVector2 = vector2 / LA.norm(vector2)

    dotProduct = np.dot(uVector1,uVector2)
    return 180 * np.arccos(dotProduct) / math.pi

def findArc(p1,p2, alpha):
    # find vector to midpoint
    pa = (p2 - p1) / 2.0
    # find rhombus center
    rhombusCenter = p1 + pa
    # get distances a and b
    a = LA.norm(pa)
    b = math.sqrt(pow(alpha,2) - pow(a,2))
    
    # distance to rhombus center
    xRC = b * (p2[1] - rhombusCenter[1]) / a
    yRC = -b * (p2[0] - rhombusCenter[0]) / a
    
    C = np.array([rhombusCenter[0] + xRC,rhombusCenter[1] + yRC])
    # find vectors from circle center
    CP1 = p1 - C
    CP2 = p2 - C
    Cx = np.array([alpha,0])
    # find arc angles
    if CP2[1] > -0.05:
        arcStartAngle = findAngle(CP2,Cx) # angle from x-axis to P2
        angleInbetween = findAngle(CP2,CP1) # angle from P2 to P1
    else:
        arcStartAngle = -findAngle(CP2,Cx) # angle from x-axis to P2
        angleInbetween = findAngle(CP2,CP1) # angle from P2 to P1
    return C, arcStartAngle, arcStartAngle + angleInbetween

########## end of functions ##########
# read YAML file
file = open("build\src\Debug\output.yaml")
parsed_yaml = yaml.load(file, Loader=yaml.FullLoader)

# Create numpy array of convex hull points
chPoints = np.array(parsed_yaml.get("convexHull_points"))

# Create numpy array of removed points and their index
removedPointsRadius = flip(np.array([i[0] for i in parsed_yaml.get("removed_points_radius_and_index")]))
removedPointsIndex = [i[1] for i in parsed_yaml.get("removed_points_radius_and_index")]
# flip the array of indexes
rIndex = np.array(removedPointsIndex[::-1])
# get array of eliminated points
removedPoints = np.array([chPoints[i] for i in removedPointsIndex])
#get the amount of eliminated points
noRemovedPoints = parsed_yaml.get("eliminated_points")

# get alpha
alpha = parsed_yaml.get("alpha")

# create the figures
fig, ax = plt.subplots()

# draw circles
n = len(chPoints)
initial_alpha = alpha

# remove the eliminated points
newSCHPoints = np.delete(chPoints,rIndex,axis=0)
n = len(newSCHPoints)

# plot points
points_plot = plt.scatter(chPoints[:,0],chPoints[:,1])
plt.scatter(removedPoints[:,0],removedPoints[:,1], color='hotpink')

# get plot limits
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()

# calculate and plot the initial arcs
for i in range(n):
    # declare the points touching the arc
    p1 = newSCHPoints[i % n]
    p2 = newSCHPoints[(i + 1) % n]

    # Arc characteristics
    C, theta1, theta2 = findArc(p1,p2,initial_alpha)
    
    # draw the arc    
    arc = mpatches.Arc(C, 2 * initial_alpha, 2 * initial_alpha, 0, theta1, theta2)
    ax.add_patch(arc)

# modify axes scale and limits
plt.subplots_adjust(bottom = 0.25)
plt.axis('scaled')
plt.xlim(xmin-1,xmax+1)
plt.ylim(ymin-1,ymax+1)

# Create horizontal slider
sliderAxes = plt.axes([0.2, 0.05, 0.6, 0.1], facecolor = 'ghostwhite')
Nslider = Slider(
    ax = sliderAxes,
    label = 'Alpha value',
    valmin = math.ceil(initial_alpha * 0.75),
    valmax = math.floor(removedPointsRadius[len(removedPointsRadius)-1]*1.25),
    valinit = initial_alpha
)

# make points plot the active axis
plt.sca(ax)

# initialize the index at 0
index = 0
# get current alpha from alpha
currentAlpha = alpha

def update(N):
    global newSCHPoints
    global n, noRemovedPoints
    global index, rIndex
    global currentAlpha
    
    # clear the axes
    plt.cla()

    # find direction of the slider
    # if true, the slider is moved from right to left,
    # else, the slider is moved form left to right
    newAlpha = N
    if newAlpha < currentAlpha:
        reverse = True
    else:
        reverse = False

    # update current alpha
    currentAlpha = newAlpha

    if reverse:
        # remove previous point from hull
        # check if alfa is smaller than the previous and current point
        if index > 0 and N < removedPointsRadius[index] and N < removedPointsRadius[index-1]:
            # delete all points from index-1 onwards
            newSCHPoints = np.delete(chPoints,rIndex[index-1:],axis=0)
          
            # update index and increase number of eliminated points
            index -= 1
            noRemovedPoints += 1
        # check if last point is active and if N is smaller than its radius
        elif noRemovedPoints == 0 and N < removedPointsRadius[index]:
            # remove last point from hull
            newSCHPoints = np.delete(chPoints,rIndex[len(rIndex)-1], axis=0)
          
            # set the amount of eliminated points to 1
            noRemovedPoints = 1
    else:
        # add points to the hull
        if N >= removedPointsRadius[index] and noRemovedPoints > 0:
            # remove all points from index+1 onwards
            newSCHPoints = np.delete(chPoints,rIndex[index+1:],axis=0)
          
            # check if index is in range
            if index < len(removedPoints) - 1:
                # if in range, increment index and decrease no. of points to add
                index += 1
                noRemovedPoints -= 1
            else:
                # else, limit to max value and set amount of points to add to 0
                index = len(removedPoints) - 1
                noRemovedPoints = 0
    
    # find new number of points in hull
    n = len(newSCHPoints)

    # calculate the arcs for all n points
    for i in range(n):
        # declare the points touched by the arc
        p1 = newSCHPoints[i % n]
        p2 = newSCHPoints[(i + 1) % n]

        # get the center, initial and ending angle
        C, theta1, theta2 = findArc(p1,p2,N)

        # create and draw the arc
        arc = mpatches.Arc(C, 2 * N, 2 * N, 0, theta1, theta2)
        ax.add_patch(arc)

    # plot points
    plt.scatter(removedPoints[:,0],removedPoints[:,1], color='hotpink')
    plt.scatter(newSCHPoints[:,0],newSCHPoints[:,1])

    # set plot limits
    plt.axis('scaled')
    plt.xlim(xmin-1,xmax+1)
    plt.ylim(ymin-1,ymax+1)

# Register the update function with each slider
Nslider.on_changed(update)
plt.show()


