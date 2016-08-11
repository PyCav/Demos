import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as anim

def run(height,width,viscosity,u0,N_t,N_rec,N_skip = 200,barrier_file = './Barriers/bow.dat'):
	# "relaxation" parameter
	omega = 1 / (3*viscosity + 0.5)
	# abbreviations for lattice-Boltzmann weight factors
	four9ths = 4.0/9.0
	one9th   = 1.0/9.0
	one36th  = 1.0/36.0
	N_anim = N_t-N_skip

	# Initialize all the arrays to steady rightward flow:
	# n 0: 0, 1: N, 2: S, 3: E, 4: W, 5: NE, 6: SE, 7: NW, 8: SW
	def initialise():
		n = np.ones((9,height,width))
		# particle densities along 9 directions
		n[0] = four9ths * (n[0] - 1.5*u0**2)

		n[1] = one9th * (n[1] - 1.5*u0**2)
		n[2] = one9th * (n[2] - 1.5*u0**2)

		n[3] = one9th * (n[3] + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
		n[4] = one9th * (n[4] - 3*u0 + 4.5*u0**2 - 1.5*u0**2)

		n[5] = one36th * (n[5] + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
		n[6] = one36th * (n[6] + 3*u0 + 4.5*u0**2 - 1.5*u0**2)

		n[7] = one36th * (n[7] - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
		n[8] = one36th * (n[8] - 3*u0 + 4.5*u0**2 - 1.5*u0**2)

		# macroscopic density
		rho = np.sum(n, axis = 0)

		# macroscopic x velocity
		ux = (n[3] + n[5] + n[6] - n[4] - n[7] - n[8]) / rho
		# macroscopic y velocity
		uy = (n[1] + n[5] + n[7] - n[2] - n[6] - n[8]) / rho

		return n, rho, ux, uy

	n, rho, ux, uy = initialise()

	# Initialize barriers:
	barrier = np.zeros((height,width), bool)

	b_map = np.loadtxt(barrier_file,dtype=bool)
	b_map_size = b_map.shape
	top_left_edge = [height/2 - b_map_size[0]/2,width/8]
	bottom_left_edge = [height/2 + b_map_size[0]/2,width/8+b_map_size[1]]
	barrier[int(top_left_edge[0]):int(bottom_left_edge[0]),int(top_left_edge[1]):int(bottom_left_edge[1])] = b_map[::-1,:]

	barrierN = np.roll(barrier,  1, axis=0)					# sites just north of barriers
	barrierS = np.roll(barrier, -1, axis=0)					# sites just south of barriers
	barrierE = np.roll(barrier,  1, axis=1)					# etc.
	barrierW = np.roll(barrier, -1, axis=1)
	barrierNE = np.roll(barrierN,  1, axis=1)
	barrierNW = np.roll(barrierN, -1, axis=1)
	barrierSE = np.roll(barrierS,  1, axis=1)
	barrierSW = np.roll(barrierS, -1, axis=1)

	# Move all particles by one step along their directions of motion (pbc):
	def stream(n):
		n[1]  = np.roll(n[1],   1, axis=0)	# axis 0 is north-south; + direction is north
		n[5] = np.roll(n[5],  1, axis=0)
		n[7] = np.roll(n[7],  1, axis=0)

		n[2]  = np.roll(n[2],  -1, axis=0)
		n[6] = np.roll(n[6], -1, axis=0)
		n[8] = np.roll(n[8], -1, axis=0)

		n[3]  = np.roll(n[3],   1, axis=1)	# axis 1 is east-west; + direction is east
		n[5] = np.roll(n[5],  1, axis=1)
		n[6] = np.roll(n[6],  1, axis=1)

		n[4]  = np.roll(n[4],  -1, axis=1)
		n[7] = np.roll(n[7], -1, axis=1)
		n[8] = np.roll(n[8], -1, axis=1)
		# Use tricky boolean arrays to handle barrier collisions (bounce-back):
		n[1][barrierN] = n[2][barrier]
		n[2][barrierS] = n[1][barrier]

		n[3][barrierE] = n[4][barrier]
		n[4][barrierW] = n[3][barrier]

		n[5][barrierNE] = n[8][barrier]
		n[8][barrierSW] = n[5][barrier]

		n[7][barrierNW] = n[6][barrier]
		n[6][barrierSE] = n[7][barrier]
		return n
		
	# Collide particles within each cell to redistribute velocities (could be optimized a little more):
	def collide(n):
		rho = np.sum(n, axis = 0)
		ux = (n[3] + n[5] + n[6] - n[4] - n[7] - n[8]) / rho
		uy = (n[1] + n[5] + n[7] - n[2] - n[6] - n[8]) / rho
		ux2 = ux * ux				# pre-compute terms used repeatedly...
		uy2 = uy * uy
		u2 = ux2 + uy2
		omu215 = 1 - 1.5*u2			# "one minus u2 times 1.5"
		uxuy = ux * uy
		n[0] = (1-omega)*n[0] + omega * four9ths * rho * omu215

		n[1] = (1-omega)*n[1] + omega * one9th * rho * (omu215 + 3*uy + 4.5*uy2)
		n[2] = (1-omega)*n[2] + omega * one9th * rho * (omu215 - 3*uy + 4.5*uy2)
		n[3] = (1-omega)*n[3] + omega * one9th * rho * (omu215 + 3*ux + 4.5*ux2)
		n[4] = (1-omega)*n[4] + omega * one9th * rho * (omu215 - 3*ux + 4.5*ux2)

		n[5] = (1-omega)*n[5] + omega * one36th * rho * (omu215 + 3*(ux+uy) + 4.5*(u2+2*uxuy))
		n[7] = (1-omega)*n[7] + omega * one36th * rho * (omu215 + 3*(-ux+uy) + 4.5*(u2-2*uxuy))
		n[6] = (1-omega)*n[6] + omega * one36th * rho * (omu215 + 3*(ux-uy) + 4.5*(u2-2*uxuy))
		n[8] = (1-omega)*n[8] + omega * one36th * rho * (omu215 + 3*(-ux-uy) + 4.5*(u2+2*uxuy))
		# Force steady rightward flow at ends (no need to set 0, N, and S components):
		n[3][:,0] = one9th * (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
		n[4][:,0] = one9th * (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)

		n[5][:,0] = one36th * (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
		n[6][:,0] = one36th * (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)

		n[7][:,0] = one36th * (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
		n[8][:,0] = one36th * (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
		return n, ux, uy, rho

	# Compute curl of the macroscopic velocity field:
	def curl(ux, uy):
		return np.roll(uy,-1,axis=1) - np.roll(uy,1,axis=1) - np.roll(ux,-1,axis=0) + np.roll(ux,1,axis=0)

	record_array = np.zeros((N_rec,height,width,2))

	j = 0 
	for i in range(N_t):
		n_stream = stream(n)
		n, ux, uy, rho = collide(n_stream)
		if i >= N_skip-1:
			if (i-N_skip-1) % int(N_anim/N_rec) == 0:
				record_array[j,:,:,0] = ux
				record_array[j,:,:,1] = ux
				j += 1

	# Here comes the graphics and animation...
	fig = plt.figure(figsize=(9,6))
	plt.tick_params(axis = 'both',which = 'both',labelbottom='off',labelleft='off')
	fluidImage = plt.imshow(ux, origin='lower', norm=plt.Normalize(-.1,.1),cmap=plt.get_cmap('jet'), interpolation='none')
	bImageArray = np.zeros((height, width, 4), np.uint8)	# an RGBA image
	bImageArray[barrier,3] = 255								# set alpha=255 only at barrier sites
	barrierImage = plt.imshow(bImageArray, origin='lower', interpolation='none')

	# Function called for each successive animation frame:
	def nextFrame(arg):							# (arg is the frame number, which we don't need)
		fluidImage.set_array(record_array[arg,:,:,0])
		return (fluidImage, barrierImage)		# return the figure elements to redraw

	animate = anim.FuncAnimation(fig, nextFrame, frames = N_rec, interval = 50, blit = True, repeat = False)
	return animate
