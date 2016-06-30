import numpy as np
import matplotlib.pyplot as plt

def decay_cycle(N_t,N,N_d,P,data):
	for i in range(1,N_t):
		rand = np.random.rand(N)
		for j in range(N):
			for k in range(N_d):
				if rand[j] < P[k]:
					if data[j,i] < N_d:
						data[j,i:] = data[j,i:] + 1.0
						break
	return data

def counter(N_t,N_d,data):
	# Counted nucleii
	count = np.zeros((N_d,N_t+1))
	for i in range(N_d):
		count[i,:] = np.sum(data.astype('int64') == i,axis = 0)
	return count

# Number of nucleii
N = 100000

# Number of daughters in chain
N_d = 4

# Probability of each decay
P = np.zeros((N_d))
P[0] = 0.25
P[1] = 0.025
P[2] = 0.075
P[3] = 0.0001

# Number of time steps
N_t = 25

# Data record array
initial_state = np.zeros((N,N_t+1))

data = decay_cycle(N_t,N,N_d,P,initial_state)

count = counter(N_t,N_d,data)

time = np.arange(0,N_t+1)
for i in range(N_d):
	plt.plot(time,count[i,:])

plt.show()
