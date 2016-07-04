from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import dblquad,quad
from shapes import *
from PIL import Image
from math import pi
from mpl_toolkits.mplot3d import Axes3D

from scipy import fftpack
import numpy as np
import pylab as py

#do not change
apertures=[]
dr=0.01

class lightSource:
    def __init__(self,color='w'):
        if(color=='r'):
            self.clr=[255,0,0]
        elif(color=='g'):
            self.clr=[0,255,0]
        elif(color=='b'):
            self.clr=[0,0,255]
        else:
            self.clr=[255,255,255]

class screen:
    def __init__(self,sizeX,sizeY,distance,color=[0,0,0]):
        self.width=sizeX
        self.height=sizeY
        self.L=distance
        self.clr=color

class aperture:
    def __init__(self,aType,r0,dim1,dim2=0.0,dim3=0.0, Rotate=0.0, invPhase=False):
        self.centre=r0
        self.dim1=dim1
        self.dim2=dim2
        self.dim3=dim3
        if(dim1>=dim2 and dim1>=dim3):
            self.extentX=dim1
            self.extentY=dim1
        elif(dim2>dim1 and dim2>=dim3):
            self.extentX=dim2
            self.extentY=dim2
        else:
            self.extentX=dim3
            self.extentY=dim3
    
        if(invPhase):
            self.phase=-1.0
        else:
            self.phase=1.0
        
        if(aType=="Circle"):
            self.aTYPE='c'
            self.shape=ellipse(r0, self.dim1,angle=Rotate)
        elif(aType=="Rectangle"):
            self.aTYPE='r'
            self.shape=rectangle(r0, self.dim1, self.dim2,angle=Rotate)
        elif(aType=="Triangle"):
            self.aTYPE='t'
        elif(aType=="Grating"):
            self.aTYPE='g'
        elif(aType=="Polygon"):
            self.aTYPE='p'
        elif(aType=="Ellipse"):
            self.aTYPE='e'
            self.shape=ellipse(r0, self.dim1, self.dim2,angle=Rotate)
        else:
            self.aTYPE='o' 
    def isInside(self,x,y):
        return self.shape.enclosesPoint(x,y)

class apertureArray:
    def __init__(self,aType,r0,maxX,maxY,spacingX,spacingY,dim1,dim2=0.0,dim3=0.0,angRotate=0.0,invPhase=False):
        if(invPhase):
            self.phase=-1.0
        else:
            self.phase=1.0
        
        if(aType=="Circle"):
            self.aTYPE='c'
        elif(aType=="Rectangle"):
            self.aTYPE='r'
        elif(aType=="Triangle"):
            self.aTYPE='t'
        elif(aType=="Grating"):
            self.aTYPE='g'
        elif(aType=="Polygon"):
            self.aTYPE='p'
        elif(aType=="Ellipse"):
            self.aTYPE='e'
        else:
            self.aTYPE='o'  
        self.apetures=[]
        self.centre=r0,
        self.extentX = maxX
        self.extentY = maxY
        
    def isInside(self,x,y):
        flag=False
        for shape in self.apertures:
            if(shape.enclosesPoint((x,y))):
                flag=True
                break
        return flag

class diffractionPattern:
    def __init__(self,dType="Fraunhofer"):
        if(dType=="Fraunhofer"):
            self.fraun=True
        else:
            self.fraun=False
        self.isH=False
    def boundary(self):
        self.boundX=[]
        self.boundY=[]
        maxX=[0.0,0.0]
        maxY=[0.0,0.0]
        for A in apertures:
            if(A.centre[0]+2.0*A.extentX>=maxX[1]):
                maxX[1]=A.centre[0]+2.0*A.extentX
            if(A.centre[0]-2.0*A.extentX<=maxX[0]):
                maxX[0]=A.centre[0]-2.0*A.extentX
            if(A.centre[1]+2.0*A.extentY>=maxY[1]):
                maxY[1]=A.centre[1]+2.0*A.extentY
            if(A.centre[1]-2.0*A.extentY<=maxY[0]):
                maxY[0]=A.centre[1]-2.0*A.extentY
        self.mx=max(abs(i) for i in maxX)
        self.my=max(abs(i) for i in maxY)
        for i in range(0,int((2*self.mx)/dr)+1):
            self.boundX.append(-abs(self.mx)+dr*i)
        for i in range(0,int((2*self.my)/dr)+1):
            self.boundY.append(-abs(self.my)+dr*i)

    def apertureFunction(self):
        self.boundary()
        flag=False
        w, d = len(self.boundX), len(self.boundY)
        self.h = [[None] * w for i in range(d)]
        for i in range(0,w):
            for j in range(0,d):
                flag=False
                for A in apertures:
                    if(A.isInside(self.boundX[i],self.boundY[j])):
                        flag=True
                if(flag):
                    try:
                        self.h[i][j]=1.0
                    except IndexError:
                        pass
                else:
                    try:
                        self.h[i][j]=0.0
                    except IndexError:
                        pass
    
    def displayAperture(self): #problems away from (0,0) with cutting off of shapes
        if(not self.isH):
            self.apertureFunction()
            isH=True
        aperturePlane = Image.new( 'RGB', (len(DF.boundX),len(DF.boundY)), "black")
        pixels = aperturePlane.load()
        for i in range(aperturePlane.size[0]):
            for j in range(aperturePlane.size[1]):
                try:
                    val=self.h[i][j]
                except IndexError:
                    val=0.0
                if(val!=0.0):
                    pixels[i,j] = (255, 255, 255)
        aperturePlane.show()
        F1 = fftpack.fft2(aperturePlane)
        F2 = fftpack.fftshift( F1 )
        psf2D = np.abs( F2 )**2
        py.figure(1)
        py.clf()
        py.imshow( psf2D )
        py.show()

    def H(self,x,y):
        i=int(round((x+self.mx)/dr))
        j=int(round((y+self.my)/dr))
        if(0<=i<len(self.h)):
            if(0<=j<len(self.h[0])):
                return H[i][j]
            else:
                return 0.0
        else:
            return 0.0

    def psi(x,y):
        pass

    def intensity(self,lightSrc):
        self.clr=lightSrc.clr
        if(not self.isH):
            self.apertureFunction()
            isH=True
        self.I=[]
        self.I = [(255.0*x / max(self.I))%255 for x in self.I]
SCREEN=screen(20,20,10,[0,0,0])
SOURCE=lightSource([255,255,255])
AP=aperture("Ellipse",[0,0],1,1,Rotate=0)
#AP2=aperture("Rectangle",[0,2],2,2,Rotate=0)
#AP3=aperture("Ellipse",[2,0],2,Rotate=0)
apertures.append((AP))
#apertures.append((AP2))
#apertures.append((AP3))
DF=diffractionPattern()
DF.intensity(SOURCE)
DF.displayAperture()




