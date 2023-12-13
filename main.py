from header import *
from func import *
#################################################################
options = {'c1': 0.5, 'c2': 0.3, 'w':0.9}
initial_guess_1 = kasigma
dimensions = initial_guess_1.size
initial_guess = initial_guess_1.reshape((1, dimensions))
n_particles=2
init_pos = np.tile(initial_guess, (n_particles, 1))
bounds=(kasigmin,kasigmax)
optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,
		dimensions=2,options=options, init_pos=init_pos, bounds=bounds)
stats = optimizer.optimize(cost_func, iters=10, expt_data=expt_data,verbose=True)
print(stats)
#################################################################