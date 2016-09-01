import matplotlib
matplotlib.use('TkAgg')

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from time import time

from matplotlib.figure import Figure

import sys
if sys.version_info[0] < 3:
	import Tkinter as Tk
else:
	import tkinter as Tk

from fluid_module import *

class Interface():

	def __init__(self,root):
		self.N_x  = 128
		self.N_y  = 128
		define_N(self.N_x,self.N_y)

		self.dens,self.dens_prev,self.u,self.v,self.u_prev,self.v_prev = initialise()

		self.diff = 0.001
		self.visc = 0.0001
		self.stoc = 50
		self.dt   = 10**-2

		self.root = root
		self.root.wm_title("Navier-Stokes Solver")

		f = Figure(figsize=(5, 4), dpi=100)
		self.ax1 = f.add_subplot(121)
		self.ax2 = f.add_subplot(122)
		self.ax1.set_xticks([])
		self.ax1.set_yticks([])
		self.ax2.set_xticks([])
		self.ax2.set_yticks([])
		self.ax2.set_aspect('equal')
		self.img = self.ax1.imshow(self.dens[1:-1,1:-1].T,cmap = 'Greys',origin='lower',vmin = 0, vmax = 2)
		self.qiv = self.ax2.quiver(self.u[1:-1,1:-1][::4,::4].T,self.v[1:-1,1:-1][::4,::4].T,pivot='mid', color='r',scale = 15)

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

		self.n_label = Tk.Label(master=self.root, text='Resolution:')
		self.n_label.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.n_entry = Tk.Entry(master=self.root, width = 8)
		self.n_entry.insert(0, str(self.N_x))
		self.n_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.e_label = Tk.Label(master=self.root, text='Diffusion:')
		self.e_label.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.e_entry = Tk.Entry(master=self.root, width = 8)
		self.e_entry.insert(0, str(self.diff))
		self.e_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.r_label = Tk.Label(master=self.root, text='Viscosity:')
		self.r_label.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.r_entry = Tk.Entry(master=self.root, width = 8)
		self.r_entry.insert(0, str(self.visc))
		self.r_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.s_label = Tk.Label(master=self.root, text='Stochastic parameter:')
		self.s_label.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.s_entry = Tk.Entry(master=self.root, width = 4)
		self.s_entry.insert(0, str(self.stoc))
		self.s_entry.pack(side=Tk.LEFT, pady = 10, padx = 5)

		self.i_button = Tk.Button(master=self.root, text='Initialise', command=self._initialise)
		self.i_button.pack(side=Tk.LEFT, padx = 20)

		Tk.mainloop()

	def _update_plot(self):
		self.img.set_data(self.dens[1:-1,1:-1].T)
		self.qiv.set_UVC(self.u[1:-1,1:-1][::4,::4].T,self.v[1:-1,1:-1][::4,::4].T)

		self.canvas.draw()

	def _initialise(self):
		N_x_prev = self.N_x
		try:
			self.N_x,self.N_y = abs(int(self.n_entry.get())),abs(int(self.n_entry.get()))
			define_N(self.N_x,self.N_y)
		except:
			self.N_x,self.N_y = self.N_x,self.N_y
		try:
			self.diff = abs(float(self.e_entry.get()))
		except:
			self.diff = 0.
		try:
			self.visc = abs(float(self.r_entry.get()))
		except:
			self.visc = 0.
		try:
			self.stoc = abs(float(self.s_entry.get()))
		except:
			self.stoc = 10.
		self.running = True
		self.dens,self.dens_prev,self.u,self.v,self.u_prev,self.v_prev = initialise()
		if N_x_prev != self.N_x:
			self.ax1.clear()
			self.ax2.clear()
			self.ax1.set_xticks([])
			self.ax1.set_yticks([])
			self.ax2.set_xticks([])
			self.ax2.set_yticks([])
			self.img = self.ax1.imshow(self.dens[1:-1,1:-1].T,cmap = 'Greys',origin='lower',vmin = 0, vmax = 2)
			self.qiv = self.ax2.quiver(self.u[1:-1,1:-1][::4,::4].T,self.v[1:-1,1:-1][::4,::4].T,pivot='mid', color='r',scale = 15)
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
			self.dens,self.u_prev,self.v_prev = forces(self.dens,self.u_prev,self.v_prev,self.stoc)
			self.u, self.v, self.u_prev, self.v_prev = vel_step(self.u, self.v, self.u_prev, self.v_prev, self.visc, self.dt)
			self.dens, self.dens_prev = dens_step(self.dens, self.dens_prev, self.u, self.v, self.diff, self.dt)
			self._update_plot()
			self.root.after(10,self._animate)

def main():
  
	root = Tk.Tk()
	root.geometry('1000x600')
	app = Interface(root)
	root.mainloop()  

if __name__ == '__main__':
	main()  