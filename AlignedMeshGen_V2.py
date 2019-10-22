import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
import csv

def DoRotation(xspan, yspan, RotRad=0):
    """Generate a meshgrid and rotate it by RotRad radians."""

    # Clockwise, 2D rotation matrix
    RotMatrix = np.array([[np.cos(RotRad),  np.sin(RotRad)],
                          [-np.sin(RotRad), np.cos(RotRad)]])

    #x, y = np.meshgrid(xspan, yspan)
    return np.einsum('ji, mni -> jmn', RotMatrix, np.dstack([x, y]))

def DetermineIntersection(x1,y1,x2,y2,x3,y3,x4,y4):
    
    t=( (x1-x3)*(y3-y4) - (y1-y3)*(x3-x4) ) / ( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) )
    u=( (x1-x2)*(y3-y4) - (y1-y2)*(x1-x3) ) / ( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) )
    
    return t,u

def reduceSize(elemLenXceil,elemLenYceil,specLenX,specLenY,Xnew,Ynew):

    toleX=elemLenXceil*2.0
    toleY=elemLenYceil*3.0
    
    searchX = toleX+specLenX
    searchY = toleY+specLenY
    
    Ynew = Ynew[Xnew >= -searchX/2]
    Xnew = Xnew[Xnew >= -searchX/2] 
    
    Ynew = Ynew[Xnew <= searchX/2]
    Xnew = Xnew[Xnew <= searchX/2] 
    
    Xnew = Xnew[Ynew >= -searchY/2] 
    Ynew = Ynew[Ynew >= -searchY/2]
    
    Xnew = Xnew[Ynew <= searchY/2] 
    Ynew = Ynew[Ynew <= searchY/2]
    
    return Xnew,Ynew

def DetermineIntersectionPoint(x1,y1,x2,y2,x3,y3,x4,y4,t,u):
    
    P1x=x1+t*(x2-x1)
    P1y=y1+t*(y2-y1)
    P2x=x3+u*(x4-x3)
    P2y=y3+u*(y4-y3)   
    
    return P1x,P1y,P2x,P2y

# Inputs
#plyStack=[45,-1,-45,-1,45,-1,-45,-45,-1,45,-1,-45,-1,45]
#plyStack=[45,-2,-1,-2,-45,-45,-2,-1,-2,45]
plyStack=[45]
elemLenX=0.1 # desired element size
fiberAngle=45*np.pi/180 # fiberAngle < 90deg
specLenX=1.5
specLenY=1.0
meshAspectRatio=0.75

nx=int(np.ceil(specLenX/elemLenX))
elemLenXceil=specLenX/nx
ny=int(np.ceil(specLenY/(elemLenXceil*meshAspectRatio)))
elemLenYceil=specLenY/ny

plt.figure(figsize=(7,7))
plt.axis([-1, 1, -1, 1])

x=np.array([-specLenX/2,specLenX/2,specLenX/2,-specLenX/2,-specLenX/2])
y=np.array([-specLenY/2,-specLenY/2,specLenY/2,specLenY/2,-specLenY/2])
plt.plot(x,y)
# -----------
specLenXp=specLenX*np.cos(fiberAngle)+specLenY*np.sin(fiberAngle)
specLenYp=specLenX*np.sin(fiberAngle)+specLenY*np.cos(fiberAngle)


nx=int(np.ceil(specLenXp/elemLenX))
elemLenXceil=specLenX/nx
ny=int(np.ceil(specLenYp/(elemLenXceil*meshAspectRatio)))
elemLenYceil=specLenY/ny

xr = np.linspace(-specLenXp/2, specLenXp/2, nx)
yr = np.linspace(-specLenYp/2, specLenYp/2, ny)

xtrial, ytrial = np.meshgrid(xr, yr)

A= DoRotation(xr, yr, -fiberAngle)

#plt.scatter(A[0],A[1])

#Xnew,Ynew = A[0].flatten(), A[1].flatten()

Xnew = np.array([0,-1])
Ynew = np.array([-1,0])

#Xnew,Ynew = reduceSize(elemLenXceil,elemLenYceil,specLenX,specLenY,Xnew,Ynew)

plt.scatter(Xnew,Ynew, s=1)

for i in range(0,len(x)-1):
    for j in range(0,len(Xnew)-1):
        x1=x[i]
        y1=y[i]        
        x2=x[i+1]
        y2=y[i+1]
        x3=Xnew[j]
        y3=Ynew[j]
        x4=Xnew[j+1]
        y4=Ynew[j+1]
        t,u = DetermineIntersection(x1,y1,x2,y2,x3,y3,x4,y4)

        if t>=0 and t<=1:  
            if u>=0 and u<=1: 
                print(t,u)
                P1x,P1y,P2x,P2y = DetermineIntersectionPoint(x1,y1,x2,y2,x3,y3,x4,y4,t,u)
                plt.scatter(P1x,P1y, s=100)
                plt.scatter(P2x,P2y, s=100)

