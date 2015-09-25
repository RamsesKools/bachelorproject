#param.py
#Parameters for the simulation

#name of this group of simulations
SimName = 'RK_1patch'

#defines the number of cycles
Ncycle1 = 1000
Ncycle2 = 1000

# defines which epsilons to simulate
eps =[]
for x in range(1,4):
	eps.append(x)

# defines which deltas to simulate
delta =[]
for x in range(1,4):
	delta.append(x)

# defines how many sites to simulate
Nsites = 10
# defines the maximum number of states in the program
Maxinterfaces = 15

filename = "lambda_test.inp"


# defines how the lambda file should be made (works for lambda_b and lambda_u)
def create_lambda_file (Epsilon, filename):
	f = open(filename, 'w')

	eps = Epsilon * 0.3
	step = (Epsilon - eps) / (Maxinterfaces)

	#if step < 1:
	#	step = 1

	Ninterfases = ((Epsilon - eps) / step)

	f.write('%f \n' % Ninterfases)
	f.write('%f \n' % eps)

	for x in range(1,int(Ninterfases)):
		f.write('%f \n' % eps)
		eps += step

	f.close()
	return

def create_lambda_u ():
	f = open("lambda_u.inp", "w")

	f.write('8 \n0 \n0 \n 0.000000001 \n0.000001 \n0.001 \n0.1 \n1.0 \n1.0\n 2.0\n')

	f.close()
	return