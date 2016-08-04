import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from mpl_toolkits.mplot3d.axes3d import Axes3D

from IPython.display import display
import ipywidgets as widgets
from ipywidgets import interact, FloatSlider

import pycav.display as dis

def rk4(R,V,a,h,ring_r,ring_steps,q):
	# Function to time step using the Runge-Kutta 4th order algorithm
	# V - velocity vector of the particle(s)
	# R - position vector of the particle(s)
	# a - function which outputs the acceleration of the particle
	# h - time step size
	# ring_r,ring_steps,q - parameters used in the a() function

	# First
	r_1 = R
	v_1 = V
	a_1 = a(r_1,ring_r,ring_steps,q)

	# Second
	r_2 = r_1 + 0.5*v_1*h
	v_2 = v_1 + 0.5*a_1*h
	a_2 = a(r_2,ring_r,ring_steps,q)

	# Third
	r_3 = r_1 + 0.5*v_2*h
	v_3 = v_1 + 0.5*a_2*h
	a_3 = a(r_3,ring_r,ring_steps,q)

	# Fourth
	r_4 = r_1 + v_3*h
	v_4 = v_1 + a_3*h
	a_4 = a(r_4,ring_r,ring_steps,q)

	# Final position and velocity
	r_f = R + (h/6.0)*(v_1+2*v_2+2*v_3+v_4)
	v_f = V + (h/6.0)*(a_1+2*a_2+2*a_3+a_4)

	return r_f,v_f


def shm_freq(ring_r,q):
	return (q/(4*np.pi*ring_r**3))**0.5

def E_field(z_lim,z_step,q,ring_r):
	# Function to plot the electric field along the z axis at x,y = 0,0
	num_p = int(z_lim*ring_r/z_step)
	z = np.linspace(-z_lim/2.0*ring_r,z_lim/2.0*ring_r,num_p)
	
	E_z = -(q/(4*np.pi))*z*(ring_r**2+z**2)**(-1.5)

	fig_E = plt.figure()
	plt.plot(z,E_z,'b-')
	plt.plot(z,-shm_freq(ring_r,q)**2*z,'r-')
	plt.xlim([np.min(z),np.max(z)])
	plt.ylim([np.min(E_z)-0.1*np.max(E_z),np.max(E_z)+0.1*np.max(E_z)])
	plt.xlabel('Distance z / (arb.units)')
	plt.ylabel('Electrostatic Force / (arb. units)')
	plt.legend(('Analytic solution','SHM approx'))
	plt.show()

class Animation(object):

	def __init__(self,a,q,ring_r,ring_steps,h,N,N_output):
		# Initialising values from the presets given at the end of this file
		self.a = a
		self.q = q
		self.ring_r = ring_r
		self.ring_steps = ring_steps
		self.h = h
		self.N = N
		self.N_output = N_output
		
		# Initialises the tracer to be off
		self.trace = False
        
		# Set up linspace for time in animation
		self.time = np.linspace(0.0,self.h*self.N,int(self.N/self.N_output))

	def create_plot(self):
		# Creates a figure which is twice as wide as it is tall
		self.fig = plt.figure(figsize = (9,4.5))

		# First subplot - 3D plot which tracks the particle motion
		self.ax1 = self.fig.add_subplot(1,2,1, projection='3d')
		# Setting the x and y axis limits
		self.ax1.set_xlim3d([-self.ring_r,self.ring_r])
		self.ax1.set_ylim3d([-self.ring_r,self.ring_r])
		# For the 2 fixed particles case, zoom in on the y axis		
		if self.ring_steps == 2:
			self.ax1.set_ylim3d([-1.5*self.R_i[2],1.5*self.R_i[2]])
		self.ax1.set_zlim3d([-1.5*self.R_i[2],1.5*self.R_i[2]])
		# Draw the ring on this subplot
		self.draw_ring()

		# Second subplot - 2D plot showing the displacements as a function of time
		self.ax2 = self.fig.add_subplot(1,2,2)
		z_max = self.R_i[2]
		self.ax2.set_xlim([0.0,self.h*self.N])
		if self.ring_steps == 2:
			# For 2 fixed particle case, orbit is not confined to within the initial particle height
			self.ax2.set_ylim([-2.0*z_max,2.0*z_max])
		else:
			self.ax2.set_ylim([-z_max,z_max])

		# Create lines to be later referenced for the particle and the tracer
		self.particle = self.ax1.plot([self.R_i[0]],[self.R_i[1]],[self.R_i[2]],'ro')[0]
		self.tracer = self.ax1.plot([self.R_i[0]],[self.R_i[1]],[self.R_i[2]],'r-.')[0]

		# Create lines to be later referenced for the 2 lines on the 2D plot
		self.z1_plot = self.ax2.plot([0.0],[z_max],'b-')[0]
		self.z2_plot = self.ax2.plot([0.0],[z_max],'r-')[0]

	def draw_ring(self):
		# Draws equally spaced in angle sections of the ring
		for i in range(self.ring_steps):
			phi_h = 2*np.pi/float(self.ring_steps)
			x = self.ring_r*np.cos(i*phi_h)
			y = self.ring_r*np.sin(i*phi_h)
			self.ax1.plot([x],[y],[0.0],'bo')

	def create_sliders(self,R_i,V_i):
		# Records the initial position and velocity
		self.R_i = R_i
		self.V_i = V_i

		self.run_button = widgets.Button(description="Run (with new values)")
		display(self.run_button)
		self.run_button.on_click(self.mainloop)

		self.trace_button = widgets.Button(description="Tracer On/Off")
		display(self.trace_button)
		self.trace_button.on_click(self.t_onoff)

		widgets.interact(self.q_update,val = FloatSlider(min=1000., max=10000., step=0.1, value=self.q, description='Charge Factor'))
		widgets.interact(self.r_update,val = FloatSlider(min=0.01, max=1., step=0.001, value=self.R_i[2], description='Particle Height'))


	def run_anim(self,R_rec,V_rec):
		# Retrieve the data produced by the RK4 routine
		self.R = R_rec
		self.V = V_rec
		# Adjusting limits in case they have been changed by slider changes
		if self.ring_steps == 2:
			self.ax1.set_ylim3d([-1.5*self.R_i[2],1.5*self.R_i[2]])
			self.ax2.set_ylim([-2.0*self.R_i[2],2.0*self.R_i[2]])
		else:
			self.ax2.set_ylim([-self.R_i[2],self.R_i[2]])
			self.shm_f = shm_freq(self.ring_r,self.q)
		self.ax1.set_zlim3d([-1.5*self.R_i[2],1.5*self.R_i[2]])

		# Animation functions
		def nextframe(arg):
			# Update particle position
			self.particle.set_data(self.R[0:2,arg])
			self.particle.set_3d_properties(self.R[2,arg])
			
			# Update tracer (if on)
			if self.trace:
				self.tracer.set_data(self.R[0:2,:arg])
				self.tracer.set_3d_properties(self.R[2,:arg])
			else:
				self.tracer.set_data(self.R_i[0:2])
				self.tracer.set_3d_properties(self.R_i[2])

			frame_num = arg+1
			
			if self.ring_steps == 2:
				# Plot z and y displacements as functions of time
				self.z1_plot.set_data(self.time[:frame_num],self.R[2,:frame_num])
				self.z2_plot.set_data(self.time[:frame_num],self.R[1,:frame_num])
			else:
				# Plot z and SHM approx as functions of time
				self.z1_plot.set_data(self.time[:frame_num],self.R[2,:frame_num])
				self.z2_plot.set_data(self.time[:frame_num],self.R[2,0]*np.cos(self.shm_f*self.time[:frame_num]))

		# Create animation object which calls nextframe function in 100ms intervals
		self.animate = anim.FuncAnimation(self.fig, nextframe, int(self.N/self.N_output), interval = 100, blit = False)
		print('Run cell below to create animation inline')

	def t_onoff(self,val):
		# On button press, change values of trace bool
		if self.trace:
			self.trace = False
			print('Tracer Off')
		else:
			self.trace = True
			print('Tracer On')

	def q_update(self,val):
		# Slider value changes
		self.q = val
        
	def r_update(self,val):
		self.R_i[2] = val

	def mainloop(self,event):
		self.create_plot()
		print('Commencing RK4 routine...')
		
		R = np.zeros((3))
		V = np.zeros((3))
		R = self.R_i
		V = self.V_i

		if self.ring_steps == 2:
			# Give particle velocity required to perform circular orbit (if SHM approx holds)
			V[1] = self.R_i[2]*shm_freq(self.ring_r,self.q)

		# Record arrays:
		R_rec = np.zeros((3,int(self.N/self.N_output)))
		V_rec = np.zeros((3,int(self.N/self.N_output)))
		i_rec = 0
		for i in range(self.N):
			R,V = rk4(R,V,self.a,self.h,self.ring_r,self.ring_steps,self.q)
			# Only record every N_output result from rk4() for efficiency
			if (i%self.N_output) == 0:
				R_rec[:,i_rec] = R
				V_rec[:,i_rec] = V
				i_rec = i_rec + 1	

		print('RK4 time stepping has been completed')

		self.run_anim(R_rec,V_rec)

