from numpy.lib.function_base import angle
from numpy.lib.scimath import arccos
import yaml
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
from matplotlib.widgets import Slider, Button
import numpy as np
from numpy import linalg as LA

########## functions ##########
def findCircumCircleRadius(a,b,c):
    A = LA.norm(b-a)
    B = LA.norm(c-b)
    C = LA.norm(a,c)
    s = (A + B + C) / 2.0
    return (A * B * C) / (4 * math.sqrt(s * (s - A) * (s - B) * (s - C)))

def findAngle(vector1, vector2):
    uVector1 = vector1 / LA.norm(vector1)
    uVector2 = vector2 / LA.norm(vector2)

    dotProduct = np.dot(uVector1,uVector2)
    return 180 * np.arccos(dotProduct) / math.pi

def findArc(p1,p2, r):
    pa = (p2 - p1) / 2.0
    rhombusCenter = p1 + pa
    a = LA.norm(pa)
    b = math.sqrt(pow(r,2) - pow(a,2))
    
    # distance to rhombus center
    xRC = b * (p2[1] - rhombusCenter[1]) / a
    yRC = -b * (p2[0] - rhombusCenter[0]) / a
    
    C = np.array([rhombusCenter[0] + xRC,rhombusCenter[1] + yRC])
    # find vectors from circle center
    CP1 = p1 - C
    CP2 = p2 - C
    Cx = np.array([r,0])
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

file = open("build\src\Debug\output.yaml")
parsed_yaml = yaml.load(file, Loader=yaml.FullLoader)

# Create numpy array of convex hull points
chPoints = np.array(parsed_yaml.get("convexHull_points"))

# create the figure
fig, ax = plt.subplots()
points_plot = plt.scatter(chPoints[:,0],chPoints[:,1])

# draw circles
n = len(chPoints)
initial_radius = 1.5
# get plot limits
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()

for i in range(n):
    p1 = chPoints[i % n]
    p2 = chPoints[(i + 1) % n]
    # Arc characteristics
    C, theta1, theta2 = findArc(p1,p2,initial_radius)
    width = 2 * initial_radius 
    height = 2 * initial_radius 
    # print("Arc " + str(i))
    # print(width,height, theta1, theta2)
    
    arc = mpatches.Arc(C, width, height, 0, theta1, theta2)
    ax.add_patch(arc)

plt.subplots_adjust(bottom = 0.25)

# Horizontal slider
sliderAxes = plt.axes([0.25, 0.05, 0.65, 0.1], facecolor = 'lightgoldenrodyellow')
Nslider = Slider(
    ax = sliderAxes,
    label = 'Number of samples',
    valmin = 1.5,
    valmax = 5.5,
    valinit = initial_radius
)

# make points plot the active axis
plt.sca(ax)

def update(N):
    plt.cla()
    for i in range(n):
        p1 = chPoints[i % n]
        p2 = chPoints[(i + 1) % n]
        C, theta1, theta2 = findArc(p1,p2,N)
        width = 2 * N
        height = 2 * N
        arc = mpatches.Arc(C, width, height, 0, theta1, theta2)
        ax.add_patch(arc)
    points_plot = plt.scatter(chPoints[:,0],chPoints[:,1])
    # set plot limits
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)


# Register the update function with each slider
Nslider.on_changed(update)
plt.show()


