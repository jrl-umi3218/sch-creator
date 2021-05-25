import numpy
from numpy.core.defchararray import add, index
from numpy.lib.function_base import angle, extract, flip
from numpy.lib.scimath import arccos
import yaml
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
from matplotlib.widgets import Slider, Button
import numpy as np
from numpy import linalg as LA

########## functions ##########
def findAngle(vector1, vector2):
    uVector1 = vector1 / LA.norm(vector1)
    uVector2 = vector2 / LA.norm(vector2)

    dotProduct = np.dot(uVector1,uVector2)
    return 180 * np.arccos(dotProduct) / math.pi

def findArc(p1,p2, alpha):
    pa = (p2 - p1) / 2.0
    rhombusCenter = p1 + pa
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
    # print(C, CP1, CP2)
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
removedPointsRadiusAndIndex = parsed_yaml.get("removed_points_radius_and_index")
removedPointsRadius = flip(np.array([i[0] for i in removedPointsRadiusAndIndex]))
removedPointsIndex = [i[1] for i in removedPointsRadiusAndIndex]
rIndex = np.array(removedPointsIndex[::-1])
eliminatedPoints = np.array([chPoints[i] for i in removedPointsIndex])
NEliminatedPoints = parsed_yaml.get("eliminated_points")

# get alpha
alpha = parsed_yaml.get("alpha")

# create the figure
fig, ax = plt.subplots()


# draw circles
n = len(chPoints)
initial_alpha = alpha

# check for removed points
newCHPoints = chPoints
if initial_alpha < removedPointsRadius[0]:
    newCHPoints = np.delete(chPoints,rIndex,axis=0)
n = len(newCHPoints)

points_plot = plt.scatter(chPoints[:,0],chPoints[:,1])
plt.scatter(eliminatedPoints[:,0],eliminatedPoints[:,1], color='hotpink')
# get plot limits
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()

for i in range(n):
    p1 = newCHPoints[i % n]
    p2 = newCHPoints[(i + 1) % n]
    # Arc characteristics
    C, theta1, theta2 = findArc(p1,p2,initial_alpha)
    width = 2 * initial_alpha 
    height = 2 * initial_alpha 
    # print("Arc " + str(i))
    # print(width,height, theta1, theta2)
    
    arc = mpatches.Arc(C, width, height, 0, theta1, theta2)
    ax.add_patch(arc)

plt.subplots_adjust(bottom = 0.25)

# Horizontal slider
sliderAxes = plt.axes([0.15, 0.05, 0.7, 0.1], facecolor = 'ghostwhite')
Nslider = Slider(
    ax = sliderAxes,
    label = 'Alpha value',
    valmin = initial_alpha * 0.75,
    valmax = removedPointsRadius[len(removedPointsRadius)-1]*1.25,
    valinit = initial_alpha
)

# make points plot the active axis
plt.sca(ax)

index = 0
reverse = False

def update(N):
    plt.cla()
    global newCHPoints
    global n, NEliminatedPoints
    global index, rIndex
    global reverse
    
    # check if eliminated point Radius is larger than current alpha
    if (N >= removedPointsRadius[index]) and NEliminatedPoints > 0:
        # insert the points to the hull, delete points after current index
        newCHPoints = np.delete(chPoints,rIndex[index+1:],axis=0)
        n = len(newCHPoints)
        if index == len(rIndex)-1:
            index = len(rIndex) - 1
            NEliminatedPoints = 0
        else:
            index = index + 1
            NEliminatedPoints = NEliminatedPoints - 1
    # check if alpha is larger than the previous eliminated point but smaller than the current one
    elif (NEliminatedPoints < len(eliminatedPoints)) and (N < removedPointsRadius[index - 1]) and (N <= removedPointsRadius[index]):
        # remove points from the hull, current index forward
        newCHPoints = np.delete(chPoints,rIndex[index-1:],axis=0)
        n = len(newCHPoints)
        # modify index
        if index > 0:
            index = index - 1 
            NEliminatedPoints = NEliminatedPoints + 1
        else:
            index = 0
            NEliminatedPoints = len(eliminatedPoints)

    for i in range(n):
        p1 = newCHPoints[i % n]
        p2 = newCHPoints[(i + 1) % n]
        C, theta1, theta2 = findArc(p1,p2,N)
        width = 2 * N
        height = 2 * N
        arc = mpatches.Arc(C, width, height, 0, theta1, theta2)
        ax.add_patch(arc)

    # plot points
    plt.scatter(eliminatedPoints[:,0],eliminatedPoints[:,1], color='hotpink')
    plt.scatter(newCHPoints[:,0],newCHPoints[:,1])

    # set plot limits
    plt.xlim(xmin-0.5,xmax)
    plt.ylim(ymin,ymax)


# Register the update function with each slider
Nslider.on_changed(update)
plt.show()


