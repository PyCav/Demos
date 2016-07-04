# coding: utf-8
from __future__ import division, print_function
from vpython import *
import numpy as np


def _element(vec,at):
	"""
	Finds element of vpython vector at a certain index, specified by at, because getting the index of a vpython vector wasn't working
	Parameters
	----------
	vec : vpython vector
		Vector you want to get an element from
	at : integer
		Component of vector that you want to find
	"""

	if at == 0:
		return vec.x
	elif at == 1:
		return vec.y
	elif at == 2:
		return vec.z
	else:
		return None

def _set_element(vec, at, to):
	"""
	Parameters
	----------
	vec : vpython vector
		Vector whose element you want to set
	at : integer
		Component of vector you want to set
	to : float
		Value you want to set the component to 
	Sets element of vpython vector at a certain index, specified by at, because getting the index of a vpython vector wasn't working
	"""
	if at == 0:
		vec.x = to
	elif at == 1:
		vec.y = to
	elif at == 2:
		vec.z = to

def _perpendicular_vector(vec):
	"""
	Gets an arbitary vector that is perpendicular to the vector given
	Parameters
	----------
	vec : vpython vector
		Vector whose perpendicular you want to find
	"""
	if _element(vec,at = 0) == 0 and _element(vec,at = 1) == 0:
		return vec.cross(vector(0,1,0))
	return vec.cross(vector(0,1,0))



class Particle(sphere):
	"""
	Class describing a Particle with a certain position and velocity.
	"""
	def __init__(self,pos = vector(0,0,0), velocity = vector(0,0,0), inverse_mass = 0.0, radius =0.0, color = color.red):
		"""
		Parameters
		----------
		pos : vpython vector
			Initial position of Particle
		velocity : vpython vector
			Initial velocity of Particle
		inverse_mass: float
			Inverse Mass of Particle (default = 0)
		radius : float
			Radius of Particle
		color: vpython color
			Color of particle
		"""
		sphere.__init__(self,pos = pos, velocity = velocity, radius = radius, make_trail = False, color = color)
		self.velocity = velocity
		self.inverse_mass = inverse_mass

	def increment_by(self,pos_increment, velocity_increment):
		'''
		Function to increment coordinates and velocity at the same time.
		Parameters
		----------
		pos_increment: vpython vector
			Increment for pos
		velocity_increment: vpython vector
			Increment for velocity
		'''
		self.pos += pos_increment
		self.velocity += velocity_increment


class Container(box):
	"""
	Class describing a cubic container for particles to be in
	"""

	@property
	def dimension(self):
		"""Property which stores the dimensions of the cube"""
		return self._dimension

	@dimension.setter
	def dimension(self,dimension):
		"""Setter so that the width, legth, and height are updated when the dimension property is set"""
		self._dimension = dimension
		self.length = dimension
		self.width = dimension
		self.height = dimension
		self.size_changed = True

	@property
	def surface_area(self):
		"""Property which describes the total surface area of the container"""
		return 6*(self.dimension)**2

	def __init__(self, dimension):
		"""
		Parameters
		----------
		dimension: float
			The dimension of the cube, i.e. the length, width, height of the cube
		"""
		box.__init__(self,pos = vector(0,0,0), length = dimension, height = dimension, width = dimension, opacity = 0.3)
		self.dimension = dimension
		self.size_changed = False

	def is_outside_container(self,particle):
		"""
		Returns None if is inside container, returns the index of the axis along which it is outside the cotainer if it is outside the container
		Parameters
		----------
		particle: Particle
			Particle which is being checked to see if inside container or not
		"""
		for i in range(0,3):
			if abs(_element(particle.pos,at=i))  > self.length/2 - particle.radius:
				return i
		return None




class Domain(object):
	"""
	Class to specify "domains" which the containers are split into for more efficient collision detection, allowing what seems like roughly O(n) collision detection.
	"""
	def __init__(self, pos, dimension):
		"""
		Parameters
		----------
		dimension: float
			The dimension of the cube, i.e. the length, width, height of the cube
		pos: vpython vector
			Position of "origin" vertex of cube	
		"""
		self.pos = pos
		self.dimension = dimension
		self.particles = []

	def contains(self,particle):
		"""
		Returns True or False, depending on whether the particle is in the domain or not.
		Parameters
		----------
		particle: Particle
			The particle that is being checked to see if is contained in the domain.
		"""
		relative_pos = particle.pos - self.pos
		for i in range(0,3):
			if _element(relative_pos,at=i) > self.dimension or _element(relative_pos,at=i)<0:
				return False
		return True


class System(object):
	"""
	Class which describes a system composed of a number of Particles in a container
	"""
	def __init__(self, average_speed, number_of_particles, particle_radius, container, inverse_mass = 1, collides = True, record_pressure = False, record_1d_velocities = False, record_speeds = False):
		"""
		Parameters
		----------
		average_speed: float
			Average speed of particles
		number_of_particles: integer
			Number of particles in this simulation
		particle_radius: float
			Radius of particles in this simulation
		container: Container
			Container within which this simulation will run
		inverse_mass: float
			Inverse mass of the particles in this simulation
		collides: Boolean
			Determines whether the particles will collide with each other or not. (If they don't collide, this simulates an ideal gas.
		record_pressure: Boolean
			Determines whether the pressure on the container will be recorded
		record_1d_velocities: Boolean
			Determines whether the system will record the velocities along the x component. It is recorded as an unsorted list.
		record_speeds: Boolean
			Determines whether the system will record the speeds of the particles. It is recorded as an unsorted list.
		"""
		#Setting variables
		self.collides = collides
		self.particles = []
		self.container = container
		self.average_speed = average_speed
		self.particle_radius = particle_radius
		self.inverse_mass = inverse_mass
		self.number_of_particles = number_of_particles
		self.record_pressure = record_pressure
		self.pressure = 0 #Pressure is set to 0 at the start
		self.steps = 0 #Number of steps taken, is set to 0 at the start
		self.pressure_history = [] #History of instantaneous values of pressure
		self.record_1d_velocities = record_1d_velocities
		self.one_d_velocities = np.zeros(number_of_particles)
		self.record_speeds = record_speeds
		self.speeds = np.zeros(number_of_particles)
		self.domains = []

		#For loop below adds number_of_particles random particles
		for i in range(0,number_of_particles):
			self.particles.append(self._random_particle(average_speed))
		
		#If collision detection between particles is on, the container is split into domains.
		if self.collides:
			self._setup_domains()
			self._assign_particles_to_domains()
		
		#Last particle makes a trail to make the motion obvious
		self.particles[0].make_trail = True
		self.particles[0].retain = 50
		self.particles[0].color = color.blue
		


	def simulate(self, dt = 0.01):
		"""
		Simulates the system going forwards by a time step dt by checking for collisions (with other particles if collides is True, and always with walls), and then steps time forwards. Also records pressure if record_pressure is True
		Parameters
		----------
		dt: float
			Determines the amount of time by which the system should step forwards
		"""
		#Collision detection between particles using domains: Only checks if two particles have collided if they are in the same domain.
		if self.collides:
			#Remake domains if the container's size has changed
			if self.container.size_changed == True:
				self._setup_domains()
				self._assign_particles_to_domains()
				self.container.size_changed = False
			#Assign particles to domains every 3 steps (Don't need to do every step, as particles are unlikely to move to different domains in 3 steps time. Needs changing depending on the velocity of the particles)
			if self.steps % 3 == 0:
				self._assign_particles_to_domains()
				
			#The actual collision detection
			for domain in self.domains:
				for index_1,particle_1 in enumerate(domain.particles):
					for index_2, particle_2 in enumerate(domain.particles):
						if (index_2 > index_1):
							collided = self._collided(particle_1,particle_2)
							if collided:
								self._collision(particle_1, particle_2)
		
		momenta_change = 0.
		#Go through all the particles
		for (index,particle) in enumerate(self.particles):
			#Check for collisions with walls
			wall_collision_index = self.container.is_outside_container(particle)
			if wall_collision_index is not None:
				_set_element(particle.velocity,wall_collision_index, to = -1*_element(particle.velocity,wall_collision_index))
				if _element(particle.pos,wall_collision_index) > 0:
					_set_element(particle.pos,wall_collision_index, to = _element(particle.pos,wall_collision_index) - particle.radius/10)
				else:
					_set_element(particle.pos,wall_collision_index, to = _element(particle.pos,wall_collision_index) + particle.radius/10)
					momenta_change += abs(_element(particle.velocity,wall_collision_index)/ particle.inverse_mass)
		# Record the pressure, but only after a certain number of steps have been taken, when the system will be in equilibrium
		if self.record_pressure and self.steps > 200:
			instantaneous_pressure = (momenta_change/dt)/self.container.surface_area
			self._update_pressure(instantaneous_pressure)
		# Record the 1d velocity
		if self.record_1d_velocities:
			self._update_1d_velocities()
		# Record the speeds
		if self.record_speeds:
			self._update_speeds()
		self._step_time(dt)

	def _collided(self,particle_1,particle_2):
		"""
		Checks if two particles have collided, returns true if they have, false if they haven't.
		Parameters
		----------
		particle_1: Particle
			First particle that is being checked for collision
		particle_2: Particle
			The particle that particle_1 is being checked for a collision with
		"""
		distance = (particle_2.pos - particle_1.pos).mag
		# If the two particles overlap, then a collision is detected:
		if distance <= particle_2.radius + particle_1.radius: 
			return True
		return False

	def _setup_domains(self):
		"""
		Sets up domains in the container. So that the algorithm functions, don't make the domain_size any less than 2r, and the step should be smaller than domain_size by at least 1 particle_radius. The actual values chosen here are just what seem to work well, do tinker with them
		"""
		self.domain_size = self.particle_radius * 17
		step = self.particle_radius * 15.5
		x1 = np.arange(- self.container.dimension/2,(1.05 * self.container.dimension)/2,step)
		x2 = np.arange(- self.container.dimension/2,(1.05 * self.container.dimension)/2,step)
		x3 = np.arange(- self.container.dimension/2,(1.05 * self.container.dimension)/2,step)
		x, y, z = np.meshgrid(x1, x2, x3)
		self._setup_domains_helper(self,x,y,z)
		self._assign_particles_to_domains()

	@np.vectorize
	def _setup_domains_helper(self,x,y,z):
		"""
		Helper function for setup_domains(), made so that @np.vectorize can be used for faster execution
		"""
		pos = vector(x, y, z)
		self.domains.append(Domain(pos, self.domain_size))

	def _assign_particles_to_domains(self):
		"""
		Assigns a list of particles that are in each domain to the list of domains.
		"""
		for domain in self.domains:
			domain.particles = []
			for particle in self.particles:
				if domain.contains(particle):
					domain.particles.append(particle)
					

	def _has_a_particle_at(self,pos):
		"""
		Checks whether a particle alread exists at a certain position.
		Parameters
		----------
		pos: vector
			Position that is being checked for the existence of a particle
		"""
		for particle in self.particles:
			distance = (particle.pos - pos).mag
			if distance <= particle.radius*2:
				return True
		return False

	def _random_particle(self,speed):
		"""
		Function to create a new random particle with a certain speed. 
		Parameters
		----------
		speed: float
			Speed of random particle to be created
		"""
		l = self.container.length - self.particle_radius * 2
		position = self._random_position_inside_system()
		there_is = self._has_a_particle_at(position)
		while there_is:
			position = vector(l*random() - l/2, l*random() - l/2,l*random() - l/2)
			there_is = self._has_a_particle_at(position)
		velocity = speed * norm(vector(1*random() - 0.5,1*random() - 0.5, 1*random() - 0.5))
		return Particle(pos = position, velocity = velocity, inverse_mass = self.inverse_mass, radius = self.particle_radius)

	def _random_position_inside_system(self):
		"""
		Finds a random position within the container.
		"""
		l = self.container.dimension - self.particle_radius
		return vector(l*random() - l/2, l*random() - l/2,l*random() - l/2)

	def _step_time(self,dt):
		"""
		Steps time forwards by dt. Only intended to be called by simulate()
		Parameters
		----------
		dt: float
			Size of time step to be taken
		"""
		for particle in self.particles:
			particle.increment_by(particle.velocity*dt, vector(0,0,0))
		self.steps += 1

	def _collision(self, particle_1, particle_2):
		"""
		Defines how to change velocities when two particles, particle_1 and particle_2 collide.
		Parameters
		----------
		particle_1: Particle
			One of the particles that has collided
		particle_2: Particle
			The other particle that has collided
		"""
		axis_1 = norm(particle_1.pos - particle_2.pos)
		particle_1.pos += axis_1 * (particle_1.pos - particle_2.pos).mag/20
		axis_2 = norm(_perpendicular_vector(axis_1))
		axis_3 = norm(axis_1.cross(axis_2))
		axes = [axis_1,axis_2,axis_3]
		particle_1_v_after = vector(0,0,0)
		particle_2_v_after = vector(0,0,0)
		for i in range (1,3):
			particle_1_v_after += particle_1.velocity.dot(axes[i])*axes[i]
			particle_2_v_after += particle_2.velocity.dot(axes[i])*axes[i]
		u1_parallel_axis_1 = particle_1.velocity.dot(axes[0])
		u2_parallel_axis_1 = particle_2.velocity.dot(axes[0])
		Z= ((u1_parallel_axis_1 * particle_2.inverse_mass) +(u2_parallel_axis_1 *particle_1.inverse_mass)) / (particle_1.inverse_mass + particle_2.inverse_mass)
		particle_1_v_after += axes[0] * (-1* u1_parallel_axis_1 + Z*2)
		particle_1.velocity = particle_1_v_after
		axis_1 = norm(particle_1.pos - particle_2.pos)
		particle_2_v_after += axes[0] * (-1* u2_parallel_axis_1 + Z*2)
		particle_2.velocity = particle_2_v_after


	def _update_pressure(self,instantaneous_pressure):
		"""
		Updates pressure of system.
		Parameters
		----------
		instantaneous_pressure: float
			Pressure at a certain time
		"""
		self.pressure_history.append(instantaneous_pressure)
		self.pressure = np.mean(self.pressure_history)


	def _update_1d_velocities(self):
		"""
		Update the list of 1d-velocities of particles
		"""
		for index, particle in enumerate(self.particles):
			self.one_d_velocities[index] = _element(particle.velocity, at = 0)

	def _update_speeds(self):
		"""
		Update the list of 3d speeds of particles
		"""
		for index, particle in enumerate(self.particles):
			self.speeds[index] = particle.velocity.mag