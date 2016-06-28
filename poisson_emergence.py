#here begins poisson_emergence.py
import matplotlib.pyplot as plt
from matplotlib import animation
import random
import copy

#find no. of counts in an interval
def counts(prob_photon, no_trials): #no of counts in an interval
    counts = 0
    for i in range(no_trials):
        x = random.uniform(0,1)
        if x < prob_photon:
            counts += 1
    return counts

#generate the poisson distribution
prob_photon = 3.0e-4 #probability of receiving a photon in interval (p --> 0)
no_trials = 10000 #number of trials per interval (n --> inf)
repeats = 1000 #number of intervals

results = [] #counts from each interval
results_sequence = [] #the list of results (for the animation)
random.seed()

for i in range(repeats):
    results.append(counts(prob_photon, no_trials))
    results_sequence.append(copy.copy(results))

#initialise the histogram
fig = plt.figure()
plt.xlim([0,8])
plt.hist(results_sequence[0])

#animate the histogram
def update_hist(frame_no, results_sequence):
    plt.cla()
    plt.hist(results_sequence[frame_no], align = 'mid', bins = [0, 1, 2, 3, 4, 5, 6, 7, 8])
    plt.xlabel("Counts per interval")
    plt.ylabel("Frequency")
    plt.title("Emergence of the Poisson Distribution")
    plt.xlim([0,8])

#create mpeg4
animation = animation.FuncAnimation(fig, update_hist, repeats, fargs=(results_sequence,), repeat = False)
animation.save('poisson_hist.mp4', writer = 'ffmpeg',  fps = 20)
