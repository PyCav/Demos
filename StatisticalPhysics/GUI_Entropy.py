import matplotlib
matplotlib.use('TkAgg')

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.patches as mpatches
import matplotlib.image as image

from time import sleep
import os 

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
		self.num_of_oscil = [100,50]
		self.init_energy = [12,5]
		self.spf = 1
		self.xlim = 24
		self.y = []
		self.contact = False
		self.running = True
		self.transys = 0

		self.create_systems(self.num_of_oscil,self.init_energy,delta)		

		self.root = root
		self.root.wm_title("Boltzmann Distribution")

		f = Figure(figsize=(5, 4), dpi=100)
		self.ax = f.add_subplot(111)

		dir_path = os.path.dirname(os.path.realpath(__file__))
		pic_path = dir_path+'\eqn_boltz.png'
		self.im = image.imread(pic_path)

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

		self.b_button = Tk.Button(master=self.root, text='Fit', command=self._fit)
		self.b_button.pack(side=Tk.RIGHT, padx = 2)

		self.t_button = Tk.Button(master=self.root, text='Contact', command=self._contact)
		self.t_button.pack(side=Tk.RIGHT, padx = 2)

		self.n_label = Tk.Label(master=self.root, text='Number of Boxes (1,2):')
		self.n_label.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.n_1_entry = Tk.Entry(master=self.root, width = 8)
		self.n_1_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)
		self.n_2_entry = Tk.Entry(master=self.root, width = 8)
		self.n_2_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.e_label = Tk.Label(master=self.root, text='Initial Energy (1,2):')
		self.e_label.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.e_1_entry = Tk.Entry(master=self.root, width = 5)
		self.e_1_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)
		self.e_2_entry = Tk.Entry(master=self.root, width = 5)
		self.e_2_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.r_label = Tk.Label(master=self.root, text='Steps per Frame:')
		self.r_label.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.r_entry = Tk.Entry(master=self.root, width = 8)
		self.r_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.i_button = Tk.Button(master=self.root, text='Initialise', command=self._initialise)
		self.i_button.pack(side=Tk.LEFT, padx = 20)

		Tk.mainloop()

	# Function to create a system of N oscillators from a given distribution
	def create_systems(self,N,energy,dist):
		self.systems = [dist(energy[0],N[0]),dist(energy[1],N[1])]

	# Functions to perform a random transfer of energy
	def rand_transfer(self,N,size = 2):
		return np.random.randint(N, size = size)

	def transfer(self,i,j,sys_num = None):
		if self.contact:
			sys_give = self.num_of_oscil.index(max(self.num_of_oscil))
			sys_recv = self.num_of_oscil.index(max(self.num_of_oscil))
			if i >= max(self.num_of_oscil):
				i -= max(self.num_of_oscil)
				sys_give = (sys_give+1)%2
			if j >= max(self.num_of_oscil):
				j -= max(self.num_of_oscil)
				sys_recv = (sys_recv+1)%2
			if self.systems[sys_give][i] != 0:
				self.systems[sys_give][i] -= 1
				self.systems[sys_recv][j] += 1
				if sys_give == 0 and sys_recv == 1:
					self.transys += 1
				elif sys_give == 1 and sys_recv == 0:
					self.transys -= 1
			self.swap = True						
		else:
			if self.systems[sys_num][i] != 0:
				self.systems[sys_num][i] -= 1
				self.systems[sys_num][j] += 1
				self.swap = True		

	def energy_swap(self):
		self.swap = False
		if self.contact:
			while not self.swap:
				rand_indices = self.rand_transfer(sum(self.num_of_oscil))
				self.transfer(rand_indices[0],rand_indices[1])
		else:
			for i in range(2):
				self.swap = False
				while not self.swap:
					rand_indices = self.rand_transfer(self.num_of_oscil[i])
					self.transfer(rand_indices[0],rand_indices[1],i)
 
	def xlimit(self):
		for sys in self.systems:
			if int(sys.max()) > self.xlim:
				self.xlim = int(sys.max())

	def ylimit(self):
		if self.n_1[0] != self.n_1.max() or self.n_2[0] != self.n_2.max():
			if self.n_1.max() < self.n_2.max():
				self.ylim = int(1.1*self.n_2.max())
			else:
				self.ylim = int(1.1*self.n_1.max())
		else:
			if self.ylim < self.n_1.max():
				self.ylim = int(1.1*self.n_1.max())
			elif self.ylim < self.n_2.max():
				self.ylim = int(1.1*self.n_2.max())

	def _update_plot(self):
		self.ax.clear()
		self.xlimit()
		self.n_1,_,_ = self.ax.hist(self.systems[0],np.arange(self.xlim),color = 'g')
		self.n_2,_,_ = self.ax.hist(self.systems[1],np.arange(self.xlim),color = 'b')
		self.ylimit()
		self.ax.set_xlabel('Quanta of Energy')
		self.ax.set_ylabel('Number of Boxes')
		self.ax.set_xlim(0.,self.xlim)
		self.ax.set_ylim(0.,self.ylim)
		if len(self.y) == 2:
			self.ax.plot(self.n_range[0],self.y[0],'r-')
			self.ax.plot(self.n_range[1],self.y[1],'r-')
			green_patch = mpatches.Patch(color='green', label=self.boxstr[0])
			blue_patch = mpatches.Patch(color='blue', label=self.boxstr[1])
			self.ax.legend(handles=[green_patch,blue_patch])
			self.ax.imshow(self.im, aspect='auto', extent=(0.79*self.xlim, 0.99*self.xlim, 0.6*self.ylim, 0.8*self.ylim))
		if self.contact:
			self.ax.set_title('Thermal Contact, Net Transfer (1 to 2): {}'.format(self.transys))
		else:
			self.ax.set_title('Thermal Isolation, Net Transfer (1 to 2): {}'.format(self.transys))
		
		self.canvas.draw()

	def _initialise(self):
		try:
			self.num_of_oscil = [abs(int(self.n_1_entry.get())),abs(int(self.n_2_entry.get()))]
		except:
			self.num_of_oscil = [100,50]
		try:
			self.init_energy = [abs(int(self.e_1_entry.get())),abs(int(self.e_2_entry.get()))]
		except:
			self.init_energy = [12,5]
		try:
			self.spf = abs(int(self.r_entry.get()))
		except:
			self.spf = 1
		self.xlim = 2*max(self.init_energy)
		self.y = []
		self.transys = 0
		self.running = True
		self.contact = False
		self.create_systems(self.num_of_oscil,self.init_energy,delta)
		self._update_plot()

	def _quit(self):
		self.root.quit()	 # stops mainloop
		self.root.destroy()  # this is necessary on Windows to prevent
					# Fatal Python Error: PyEval_RestoreThread: NULL tstate

	def _fit(self):
		no_fit = False
		n_list = [self.n_1,self.n_2]
		self.n_range = []
		self.y = []
		self.boxstr = []
		temp_marker = 1
		for n in n_list:

			if len(np.nonzero(n)[0]) == len(n):
				upper_bound = len(n)
			elif n[0] == 0:
				no_fit = True
			else:
				upper_bound = np.where(n == 0)[0][0]

			if not no_fit:
				n_bins = np.arange(upper_bound)+0.5
				p = np.polyfit(n_bins,np.log(n[:upper_bound]),1)
				n_bins = np.arange(self.xlim)
				self.n_range.append(n_bins)
				self.y.append(np.exp(p[0]*n_bins+p[1]))
				self.boxstr.append('T{0} = {1:3.1f}'.format(temp_marker,-1/p[0]))
				temp_marker += 1
		self._update_plot()

	def _contact(self):
		if self.contact:
			self.contact = False
		else:
			self.contact = True

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