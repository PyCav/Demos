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
		
	def random_position_inside(self):
		l = self.dimension
		return vector(l*random() - l/2, l*random() - l/2,l*random() - l/2)
		

	
class System(object):
	"""
	Class which describes a system composed of a number of Particles in a container
	"""
	def __init__(self, average_speed, number_of_particles, particle_radius, container, inverse_mass = 1, collides = True, record_pressure = False):
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
		
		#For loop below adds number_of_particles random particles
		for i in range(0,number_of_particles):
			self.particles.append(self._random_particle(average_speed))
			
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
		momenta_change = 0.
		
		#Go through all the particles
		for (index,particle) in enumerate(self.particles):
			
			#Check for collisions with other particles if collides is True
			if self.collides:
				other_particle = self._check_collisions_for_particle(index,particle)
				if other_particle is not None:
					self._collision(particle,other_particle)
					
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
		self._step_time(dt)

	
	
	def _check_collisions_for_particle(self,index,particle):
		"""
		Checks whether a particle at a certain index is colliding with any other particle. Returns None if didn't collide with anything, returns particle it collided with if it did collide
		Parameters
		----------
		index: integer
			Index of particle being checked. This is used so that this function only checks forwards from the particle, as else each pair of particles would check for collisions twice. Set to zero to check with all particles in the system
		particle: Particle
			Particle which is being checked for collisions
		"""
		for index2, other_particle in enumerate(self.particles):
			if index2 > index:
				distance = (other_particle.pos - particle.pos).mag
				# If the two particles overlap, then a collision is detected:
				if distance <= other_particle.radius + particle.radius: 
					return other_particle
		return None
		
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
		position = self.container.random_position_inside()
		there_is = self._has_a_particle_at(position)
		while there_is:
			position = vector(l*random() - l/2, l*random() - l/2,l*random() - l/2)
			there_is = self._has_a_particle_at(position)
		velocity = speed * vector(1*random() - 0.5,1*random() - 0.5, 1*random() - 0.5)
		return Particle(pos = position, velocity = velocity, inverse_mass = self.inverse_mass, radius = self.particle_radius)
	
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
			
		