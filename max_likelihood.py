from hmm import *
from nanoporereads import *
from scipy.stats import bernoulli, binom


def neg_loglik(thetas, n, xs, zs):
    return -np.sum([binom(n, thetas[z]).logpmf(x) for (x, z) in zip(xs, zs)])

"""
m = 10
theta_A = 0.8
theta_B = 0.3
theta_0 = [theta_A, theta_B]

coin_A = bernoulli(theta_A)
coin_B = bernoulli(theta_B)

xs = map(sum, [coin_A.rvs(m), coin_A.rvs(m), coin_B.rvs(m), coin_A.rvs(m), coin_B.rvs(m)])
zs = [0, 0, 1, 0, 1]


xs = np.array(xs)


ml_A = np.sum(xs[[0,1,3]])/(3.0*m)
ml_B = np.sum(xs[[2,4]])/(2.0*m)


bnds = [(0,1), (0,1)]
minimize(neg_loglik, [0.5, 0.5], args=(m, xs, zs), bounds=bnds, method='tnc', options={'maxiter': 100})
"""

thetas = self.emission
n = 322
xs = 0
zs = [snps]
z = minimize(neg_loglik(), HmmHaplotypes.emission, args=(m, xs, zs), method='nelder-mead', options={'xtol': 1e-8, 'disp': True})