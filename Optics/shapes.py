from math import sin,cos
def rotateX(ang,r,O):
	return (r[0]-O[0])*cos(ang)+(r[1]-O[1])*sin(ang)+O[0]
def rotateY(ang,r,O):
	return -(r[0]-O[0])*sin(ang)+(r[1]-O[1])*cos(ang)+O[1]

class ellipse:
	def __init__(self,r0,dim1,dim2=0,angle=0):
		self.centre=r0
		self.dimH=dim1
		if(dim2==0):
			self.dimV=dim1
		else:
			self.dimV=dim2
		self.ang=angle
	def enclosesPoint(self,x,y):
		if(self.ang!=0):
			X=rotateX(self.ang,[x,y],self.centre)
			Y=rotateY(self.ang,[x,y],self.centre)
		else:
			X=x
			Y=y
		if(round(((X-self.centre[0])**2)/(self.dimH**2)+((Y-self.centre[1])**2)/(self.dimV**2),4)<=1.0):
			return True
		else:
			return False

class rectangle:
	def __init__(self,r0,dim1,dim2,angle=0):
		self.centre=r0
		self.dimH=dim1
		self.dimV=dim2
		self.ang=angle

	def enclosesPoint(self,x,y):
		if(self.ang!=0):
			X=rotateX(self.ang,[x,y],self.centre)
			Y=rotateY(self.ang,[x,y],self.centre)
		else:
			X=x
			Y=y
		if(self.centre[0]-self.dimH/2.0<=X<=self.centre[0]+self.dimH/2.0):
			if(self.centre[1]-self.dimV/2.0<=Y<=self.centre[1]+self.dimV/2.0):
				return True
			else:
				return False
		else:
			return False
class triangle:
	def __init__(self,dim1,dim2,dim3):
		self.centre=r0
		self.dimH=dim1
		self.dimV=dim2
		self.ang=angle
	def enclosesPoint(self):
		pass
class polygon:
	def __init__(self):
		pass
	def enclosesPoint(self):
		pass
class grating:
	def __init__(self):
		pass
	def enclosesPoint(self):
		pass