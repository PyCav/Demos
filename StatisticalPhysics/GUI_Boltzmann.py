import matplotlib
matplotlib.use('TkAgg')

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from time import sleep

from matplotlib.figure import Figure

import sys
if sys.version_info[0] < 3:
	import Tkinter as Tk
else:
	import tkinter as Tk

# Set of initial energy distribution functions
def delta(energy,N):
	return energy*np.ones(N)

def uniform(energy,N):
	return np.random.randint(energy,size = N)

class Interface():

	def __init__(self,root):
		self.num_of_oscil = 100
		self.init_energy = 12
		self.spf = 1
		self.y = None
		self.xlim = 24

		self.create_system(self.num_of_oscil,self.init_energy,delta)		

		self.root = root
		self.root.wm_title("Boltzmann Distribution")

		f = Figure(figsize=(5, 4), dpi=100)
		self.ax = f.add_subplot(121)
		self.ax_grid = f.add_subplot(122)

		self.canvas = FigureCanvasTkAgg(f, master=root)
		self.canvas.show()
		self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

		self._update_plot()

		self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

		self.q_button = Tk.Button(master=self.root, text='Quit', command=self._quit)
		self.q_button.pack(side=Tk.RIGHT, padx = 2)

		self.s_button = Tk.Button(master=self.root, text='Step', command=self._step)
		self.s_button.pack(side=Tk.RIGHT, padx = 2)

		self.r_button = Tk.Button(master=self.root, text='Run/Stop', command=self._run)
		self.r_button.pack(side=Tk.RIGHT, padx = 2)
		self.running = True

		self.b_button = Tk.Button(master=self.root, text='Fit', command=self._fit)
		self.b_button.pack(side=Tk.RIGHT, padx = 2)		

		self.n_label = Tk.Label(master=self.root, text='Number of Boxes:')
		self.n_label.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.n_entry = Tk.Entry(master=self.root)
		self.n_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.e_label = Tk.Label(master=self.root, text='Initial Energy:')
		self.e_label.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.e_entry = Tk.Entry(master=self.root)
		self.e_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.r_label = Tk.Label(master=self.root, text='Steps per Frame:')
		self.r_label.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.r_entry = Tk.Entry(master=self.root)
		self.r_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.i_button = Tk.Button(master=self.root, text='Initialise', command=self._initialise)
		self.i_button.pack(side=Tk.LEFT, padx = 20)

		Tk.mainloop()

	# Function to create a system of N oscillators from a given distribution
	def create_system(self,N,energy,dist):
		self.system = dist(energy,N)

	# Functions to perform a random transfer of energy
	def rand_transfer(self,N):
		return np.random.randint(N, size = 2)

	def energy_swap(self):
		swap = False
		while not swap:
			rand_indices = self.rand_transfer(self.num_of_oscil)
			if self.system[rand_indices[0]] != 0:
				self.system[rand_indices[0]] -= 1
				self.system[rand_indices[1]] += 1
				swap = True
 
	def xlimit(self):
		if int(self.system.max()) > self.xlim:
			self.xlim = int(self.system.max())

	def ylimit(self):
		if self.n[0] != self.n.max():
			self.ylim = int(1.1*self.n.max())
		else:
			if self.ylim < self.n.max():
				self.ylim = int(1.1*self.n.max())

	def _update_plot(self):
		self.ax.clear()
		self.ax_grid.clear()
		self.ax_grid.set_xticks([])
		self.ax_grid.set_yticks([])
		self.xlimit()
		self.n,_,_ = self.ax.hist(self.system,np.arange(self.xlim),color = 'g')
		self.ylimit()
		self.ax.plot(self.system[0]+0.5,0.,'ro')
		self.ax.set_xlabel('Quanta of Energy')
		self.ax.set_ylabel('Number of Boxes')
		self.ax.set_xlim(0.,self.xlim)
		self.ax.set_ylim(0.,self.ylim)
		if self.y != None:
			self.ax.plot(self.n_range,self.y,'r-')

		if len(self.system) in [4,9,16,25,36,49,64,81,100,121,144,169,196,225,256]:

			dim = int(len(self.system)**0.5)
			data = self.system.reshape(dim,dim)
			self.ax_grid.imshow(data,cmap = 'viridis',interpolation = 'none',vmin = 0,vmax = 3*self.init_energy)
			for (i, j), z in np.ndenumerate(data):
				self.ax_grid.text(j, i, '{:d}'.format(int(z)), ha='center', va='center', fontsize = 10)
			self.ax_grid.set_xlim(-0.5,dim-0.5)
			self.ax_grid.set_ylim(-0.5,dim-0.5)

		else:
			self.ax_grid.text(0,0,'Input either non-square or too large', ha='center', va='center')
			self.ax_grid.set_xlim(-1.5,1.5)
			self.ax_grid.set_ylim(-1.5,1.5)

		self.canvas.draw()

	def _fit(self):
		no_fit = False
		if len(np.nonzero(self.n)[0]) == len(self.n):
			upper_bound = len(self.n)
		elif self.n[0] == 0:
			no_fit = True
		else:
			upper_bound = np.where(self.n == 0)[0][0]

		if not no_fit:
			self.n_range = np.arange(upper_bound)+0.5
			p = np.polyfit(self.n_range,np.log(self.n[:upper_bound]),1)
			self.n_range = np.arange(self.xlim)
			self.y = np.exp(p[0]*self.n_range+p[1])
		self._update_plot()

	def _initialise(self):
		try:
			self.num_of_oscil = abs(int(self.n_entry.get()))
		except:
			self.num_of_oscil = 100
		try:
			self.init_energy = abs(int(self.e_entry.get()))
		except:
			self.init_energy = 12
		try:
			self.spf = abs(int(self.r_entry.get()))
		except:
			self.spf = 1
		self.xlim = 2*self.init_energy
		self.y = None
		self.running = True
		self.create_system(self.num_of_oscil,self.init_energy,delta)
		self._update_plot()

	def _quit(self):
		self.root.quit()	 # stops mainloop
		self.root.destroy()  # this is necessary on Windows to prevent
					# Fatal Python Error: PyEval_RestoreThread: NULL tstate

	def _step(self):
		self.energy_swap()
		self._update_plot()

	def _run(self):
		if not self.running:
			self.running = True
		else:
			self.running = False
		self._animate()

	def _animate(self):
		if not self.running:
			for i in range(self.spf):
				self.energy_swap()
			self._update_plot()
			self.root.after(10,self._animate)

def main():
  
	root = Tk.Tk()
	root.geometry('1000x600')
	app = Interface(root)
	root.mainloop()  

if __name__ == '__main__':
	main()  