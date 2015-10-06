import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties
import numpy as np
import sys
sys.dont_write_bytecode=True
import os
import tarfile

#delete after use
import pprint

# reads the rate from rate_matrix.dat and returns it.
# input: From, To defines which rate(ex: 1, 2 gives rate from state 1 to 2)
def read_rate(From, To):
	rate_matrix = [[0 for x in range(Nstates)] for x in range(Nstates)]

	with open('rate_matrix.dat', 'r') as f:
		lines = f.readlines()
		count1 = 0
		for line in lines:
			elements = line.split()
			count2 = 0
			for element in elements:
				rate_matrix[count1][count2] = float(element)
				count2 += 1

			count1 += 1
	return rate_matrix[From-1][To-1]

# function reads jobid file and uses SimName to m
def move_to_target(SimName):
	with open('jobid', 'r') as f:
		line = f.readline()
	jobid = line.split()[-1]
	#print 'This is the jobid of the file ' + str(jobid)

	DIR = SimName + '-' + jobid
	os.chdir(DIR)

	# this grabs the rate_matrix from the tarfile
	tarf = tarfile.open("data.tar", "r")
	tarf.extract("rate_matrix.dat")
	tarf.close()
	return

def plot(X, Y, Graphs, Type, State_State):
	colors = iter(cm.rainbow(np.linspace(0, 1, len(Graphs))))
	if Type == 1:
		XTitle = 'delta'

	elif Type == 2:
		XTitle = 'epsilon'
		Y = Y.T
	else:
		print "Incorrect plot type, correct plottypes are 1 for row by columns and 2 for columns by rows."
		exit()

	Title = "Plot of rate " + str(State_State) + " for different values of " + XTitle

	#not sure if this works
	plt.figure(Title)

	for x in range(len(Graphs)):
		plt.plot(X, Y[x], label= str(Graphs[x]), color=next(colors))
	plt.ylabel('rate')
	plt.xlabel(XTitle)
	plt.grid(True)
	plt.title(Title)
	fontP = FontProperties()
	fontP.set_size('small')
	plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0., prop = fontP)
	FileName = 'rate_' + str(State_State) + '_' + XTitle +'.png'
	plt.savefig(FileName)

print 'Lets start plotting some rates!'

# the simulation name that is used in the cycletis.py and defined in param.py
SimName = 'RK_1patch'

# defines epsilons and deltas that need to be plotted
eps_range = range(10,21,1)
delta_range = range(10,65,5)
Nstates = 2
print 'Reading rate matrices'

# this creates two arrays that are used to store all the rates
# These arrays have to be split later
rate_ub = np.array([])
rate_bu = np.array([])

# reads all the rates for all eps and deltas
for eps in eps_range:
	epsdir = 'eps' + str(eps)
	if (not os.path.isdir(epsdir)):
		print "directory " + epsdir + " does not exist"
		exit()
	os.chdir(epsdir)

	#print os.getcwd()

	for delta in delta_range:
		deltadir = 'delta' + str(delta)
		if (not os.path.isdir(deltadir)):
			print "directory " + deltadir + " does not exist"
			exit()

		os.chdir(deltadir)
		#print os.getcwd()
		move_to_target(SimName)

		rate_bu = np.append(rate_bu, read_rate(1,2))
		rate_ub = np.append(rate_ub, read_rate(2,1))
		os.chdir("../../")
	os.chdir("../")

print 'Done reading rate matrices'

print 'Lets start plotting some rates'
# first split the rate lists
rate_ub_s = rate_ub.reshape(len(eps_range), len(delta_range))
rate_bu_s = rate_bu.reshape(len(eps_range), len(delta_range))

#pprint.pprint(list(rate_ub_s))
#pprint.pprint(list(delta_range))
plot(eps_range, rate_ub_s, delta_range, 1, 'u->b')
plot(eps_range, rate_bu_s, delta_range, 1, 'b->u')
plot(delta_range, rate_ub_s, eps_range, 2, 'u->b')
plot(delta_range, rate_bu_s, eps_range, 2, 'b->u')
#plot rate as as function of delta
#colors = iter(cm.rainbow(np.linspace(0, 1, len(eps_range))))
#plt.figure(0)
#for x in range(len(eps_range)):
#	plt.plot(delta_range, rate_bu_s[x], label= str(eps_range[x]), color=next(colors))
#plt.ylabel('rate u->b')
#plt.xlabel('Delta')
#plt.grid(True)
#plt.title('Plot of rate unbound->bound for different epsilons')
#fontP = FontProperties()
#fontP.set_size('small')
#plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0., prop = fontP, title='epsilons')
#plt.savefig('rate_ub_delta.png')
#
##plot rate as function of epsilon
#colors = iter(cm.rainbow(np.linspace(0, 1, len(delta_range))))
#plt.figure(1)
#for x in range(len(delta_range)):
#	plt.plot(eps_range, rate_bu_s.T[x], label= str(delta_range[x]), color=next(colors))
#plt.ylabel('rate u->b')
#plt.xlabel('Epsilon')
#plt.grid(True)
#plt.title('Plot of rate unbound->bound for different deltas')
#fontP = FontProperties()
#fontP.set_size('small')
#plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0., prop = fontP, title='Deltas')
#plt.savefig('rate_ub_epsilon.png')



# plots the different rates. X, Y and Graphs must all be lists
# Y is a numpy 2d-array containing all the rates, X is a list with all the x values
# Graphs must be a list that contains all the different corresponding column(or row) values of each line
# plot per row or per column depending on type type =1 plot per row, type = 2 plot per column
# State_State must be a string that indicates from which state to which state the rate is
