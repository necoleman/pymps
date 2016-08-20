import numpy as np
from scipy.special import jv
from scipy import linalg
import sys

def evaluate_basis_func(angle, eig, rvec, tvec, kvec, mode='Dirichlet'):
	"""Take in index vector k, vector of rs, and vector of thetas.
	Return matrix A."""
	if mode == 'Dirichlet':
		return jv(angle*kvec,np.sqrt(eig)*rvec)*np.sin(angle*kvec*tvec)
	elif mode == 'Neumann':
		return jvp(angle*kvec,np.sqrt(eig)*rvec,1)*np.cos(angle*kvec*tvec)
	else:
		print 'Sorry, that is not a valid mode. Please choose Dirichlet or Neumann.'
		return None

def find_sing_val(A, num_bdy_pts):
	"""Find minimal singular value of boundary part of the Q of A's QR."""
	Q,R = linalg.qr(A, mode='economic')
	print Q[:num_bdy_pts-1,:].shape
	return np.min(linalg.svd(Q[:num_bdy_pts,:])[1])

def compute_square_eigenvalues():
	"""Compute Laplace Dirichlet eigenvalues of the square"""
	a = 2
	N = 25
	num_pts = 2*N
	k = np.arange(1,N).reshape(N-1,1)

	# expand about origin
	x1 = np.ones(num_pts)
	y1 = np.linspace(0,1,num_pts)
	x2 = np.linspace(0,1,num_pts)
	y2 = np.ones(num_pts)

	xint = np.random.random(2*num_pts)
	yint = np.random.random(2*num_pts)

	x = np.concatenate( (x1, x2, xint) )
	y = np.concatenate( (y1, y2, yint) )

	r = np.sqrt(x**2 + y**2)
	t = np.arctan(y/x)

	lams = np.arange(0.5, 100, 0.5)
	S = []
	for lam in lams:
		A = evaluate_basis_func(a, lam, r, t, k)
		A = A.T
		S.append(find_sing_val(A,2*num_pts))

	with open('shape.dat', 'w+') as f:
		for n in range(4*num_pts):
			f.write( str(r[n]*np.cos(t[n])) + ' ' + str(r[n]*np.sin(t[n])) + '\n')

	with open('out.dat', 'w+') as f:
		for n in range(len(lams)):
			f.write( str(lams[n]) + ' ' + str(S[n]) + '\n')

	mins = []
	for i in range(len(lams))[1:-1]:
		if S[i-1] > S[i] and S[i+1] > S[i]:
			mins.append(lams[i])

	print mins

	eigs = []
	for m in range(1,5):
		for n in range(1,5):
			eigs.append( np.pi**2*(m**2 + n**2) )
	print sorted(eigs)

if __name__ == '__main__':
	compute_square_eigenvalues()
