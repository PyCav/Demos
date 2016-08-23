import matplotlib
matplotlib.use('TkAgg')

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.image as image

from scipy.optimize import curve_fit

from time import sleep

from matplotlib.figure import Figure

import os
import sys
if sys.version_info[0] < 3:
	import Tkinter as Tk
else:
	import tkinter as Tk

def uniform(low_p,high_p,N):
	return (high_p-low_p)*np.random.rand(N)+low_p

def gaussian(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

class Interface():

	def __init__(self,root):
		self.num_of_oscil = 500
		self.m = [1,1]
		self.init_max_mom = 10
		self.spf = 10
		self.y = None
		self.xlim = 24

		self.create_system(self.num_of_oscil,self.init_max_mom,uniform)		

		self.root = root
		self.root.wm_title("Maxwell Distribution")

		f = Figure(figsize=(5, 4), dpi=100)
		self.ax1, self.ax2 = f.add_subplot(121), f.add_subplot(122)

		dir_path = os.path.dirname(os.path.realpath(__file__))
		pic_path = dir_path+'\eqn_maxwell_1.png'
		self.im_1 = image.imread(pic_path)
		pic_path = dir_path+'\eqn_maxwell_2.png'
		self.im_2 = image.imread(pic_path)

		self.canvas = FigureCanvasTkAgg(f, master=root)
		self.canvas.show()
		self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

		self._update_plot()

		self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

		self.q_button = Tk.Button(master=self.root, text='Quit', command=self._quit)
		self.q_button.pack(side=Tk.RIGHT, padx = 2)

		self.r_button = Tk.Button(master=self.root, text='Run/Stop', command=self._run)
		self.r_button.pack(side=Tk.RIGHT, padx = 2)
		self.running = True

		self.b_button = Tk.Button(master=self.root, text='Fit', command=self._fit)
		self.b_button.pack(side=Tk.RIGHT, padx = 2)		

		self.n_label = Tk.Label(master=self.root, text='Number of Boxes:')
		self.n_label.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.n_entry = Tk.Entry(master=self.root, width = 8)
		self.n_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.e_label = Tk.Label(master=self.root, text='Max. Momentum:')
		self.e_label.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.e_entry = Tk.Entry(master=self.root, width = 5)
		self.e_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.m_label = Tk.Label(master=self.root, text='Masses (1,2):')
		self.m_label.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.m1_entry = Tk.Entry(master=self.root, width = 5)
		self.m1_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.m2_entry = Tk.Entry(master=self.root, width = 5)
		self.m2_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.r_label = Tk.Label(master=self.root, text='Steps per Frame:')
		self.r_label.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.r_entry = Tk.Entry(master=self.root, width = 8)
		self.r_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.i_button = Tk.Button(master=self.root, text='Initialise', command=self._initialise)
		self.i_button.pack(side=Tk.LEFT, padx = 20)

		Tk.mainloop()

	# Function to create a system of N oscillators from a given distribution
	def create_system(self,N,max_mom,dist):
		self.momenta = [dist(-max_mom,max_mom,2*N),dist(-max_mom,max_mom,2*N)]

	# Functions to perform a random transfer of energy
	def rand_transfer(self,N):
		return np.random.randint(N, size = 2)

	def rand_angle(self):
		return 2*np.pi*np.random.rand()

	def collision(self):
		i,j = self.rand_transfer(2*self.num_of_oscil)
		theta = self.rand_angle()

		if i < self.num_of_oscil:
			mi = self.m[0]
		else:
			mi = self.m[1]

		if j < self.num_of_oscil:
			mj = self.m[0]
		else:
			mj = self.m[1]

		new_p1_x = ((mj*self.momenta[0][i]-mi*self.momenta[0][j])*np.cos(theta)+
					(mj*self.momenta[1][i]-mi*self.momenta[1][j])*np.sin(theta)+
					(mi*self.momenta[0][i]+mi*self.momenta[0][j]))
		new_p1_y = ((mi*self.momenta[0][j]-mj*self.momenta[0][i])*np.sin(theta)+
					(mj*self.momenta[1][i]-mi*self.momenta[1][j])*np.cos(theta)+
					(mi*self.momenta[1][i]+mi*self.momenta[1][j]))
		new_p2_x = ((mi*self.momenta[0][j]-mj*self.momenta[0][i])*np.cos(theta)+
					(mi*self.momenta[1][j]-mj*self.momenta[1][i])*np.sin(theta)+
					(mj*self.momenta[0][i]+mj*self.momenta[0][j])) 
		new_p2_y = ((mj*self.momenta[0][i]-mi*self.momenta[0][j])*np.sin(theta)+
					(mi*self.momenta[1][j]-mj*self.momenta[1][i])*np.cos(theta)+
					(mj*self.momenta[1][i]+mj*self.momenta[0][j]))

		self.momenta[0][i] = new_p1_x/(mi+mj)
		self.momenta[1][i] = new_p1_y/(mi+mj)
		self.momenta[0][j] = new_p2_x/(mi+mj)
		self.momenta[1][j] = new_p2_y/(mi+mj)
 
	def xlimit(self):
		if self.momenta[0].max() > self.xlim:
			self.xlim = self.momenta[0].max()

	def ylimit(self):
		if self.n1.max() > self.n2.max():
			self.ylim = int(1.1*self.n1.max())
		else:
			self.ylim = int(1.1*self.n2.max())

	def _update_plot(self):
		self.ax1.clear()
		self.ax2.clear()
		self.xlimit()
		self.bins = 0.5*np.arange(-2*self.xlim,2*self.xlim)
		self.n1,_,_ = self.ax1.hist(self.momenta[0][:self.num_of_oscil],self.bins,color = 'r')
		self.n2,_,_ = self.ax2.hist(self.momenta[0][self.num_of_oscil:],self.bins,color = 'r')
		self.ylimit()
		self.ax1.set_xlabel('Number of Quanta')
		self.ax1.set_ylabel('Number of Boxes')
		self.ax1.set_xlim(-self.xlim,self.xlim)
		self.ax1.set_ylim(0.,self.ylim)
		self.ax2.set_xlabel('Number of Quanta')
		self.ax2.set_xlim(-self.xlim,self.xlim)
		self.ax2.set_ylim(0.,self.ylim)
		self.ax2.set_yticklabels([])
		self.ax1.yaxis.set_ticks_position('left')
		self.ax2.yaxis.set_ticks_position('left')
		self.ax1.xaxis.set_ticks_position('bottom')
		self.ax2.xaxis.set_ticks_position('bottom')
		if self.y != None:
			self.ax1.plot(self.bin_centres,self.y[0],'k-')
			self.ax2.plot(self.bin_centres,self.y[1],'k-')
			self.ax1.imshow(self.im_1, aspect='auto', extent=(0., 0.99*self.xlim, 0.89*self.ylim, 0.99*self.ylim))
			self.ax2.imshow(self.im_2, aspect='auto', extent=(0., 0.99*self.xlim, 0.89*self.ylim, 0.99*self.ylim))
			self.ax1.set_title('FWHM = {0:3.2f}'.format(2.355*self.sigma[0]))
			self.ax2.set_title('FWHM = {0:3.2f}'.format(2.355*self.sigma[1]))

		self.canvas.draw()

	def _fit(self):
		p0 = [self.ylim, 0., 1.]
		self.bin_centres = (self.bins[:-1] + self.bins[1:])/2
		coeff1,_ = curve_fit(gaussian, self.bin_centres, self.n1, p0=p0)
		coeff2,_ = curve_fit(gaussian, self.bin_centres, self.n2, p0=p0)
		self.y = [gaussian(self.bin_centres,*coeff1),gaussian(self.bin_centres,*coeff2)]
		self.sigma = [coeff1[2],coeff2[2]]
		self._update_plot()

	def _initialise(self):
		try:
			self.num_of_oscil = abs(int(self.n_entry.get()))
		except:
			self.num_of_oscil = 500
		try:
			self.init_max_mom = abs(float(self.e_entry.get()))
		except:
			self.init_max_mom = 10
		try:
			self.spf = abs(int(self.r_entry.get()))
		except:
			self.spf = 10
		try:
			self.m = [abs(float(self.m1_entry.get())),abs(float(self.m2_entry.get()))]
		except:
			self.m = [1,1]
		self.xlim = 2*self.init_max_mom
		self.y = None
		self.running = True
		self.create_system(self.num_of_oscil,self.init_max_mom,uniform)	
		self._update_plot()

	def _quit(self):
		self.root.quit()	 # stops mainloop
		self.root.destroy()  # this is necessary on Windows to prevent
					# Fatal Python Error: PyEval_RestoreThread: NULL tstate

	def _run(self):
		if not self.running:
			self.running = True
		else:
			self.running = False
		self._animate()

	def _animate(self):
		if not self.running:
			for i in range(self.spf):
				self.collision()
			self._update_plot()
			self.root.after(10,self._animate)

def main():
  
	root = Tk.Tk()
	root.geometry('1000x600')
	app = Interface(root)
	root.mainloop()  

if __name__ == '__main__':
	main()  