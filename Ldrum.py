import numpy as np
import scipy.special as sp
import scipy.linalg as linalg

# Compute the first 3 eigenvalues by the method of particular solutions
# (Revived) using interior and boundary points. Betcke and Trefethen, 2003
# Translated from matlab to python by Neal Coleman, 2016

# Compute subspace angles for values of lambda
N = 36
num_pts = 2*N
k = np.arange(1,N).reshape(N-1,1)
t1 = 1.5*np.pi*np.linspace(0.5,num_pts-0.5,num=num_pts)/num_pts
r1 = 1./np.maximum(np.abs(np.sin(t1)), np.abs(np.cos(t1)))

#print len(t1)
#print len(r1)

t2 = 1.5*np.pi*np.random.random(num_pts)
r2 = np.random.random(num_pts)/np.maximum(np.abs(np.sin(t2)),np.abs(np.cos(t2)))

#print len(t2)
#print len(r2)

#print k
#print t1
#print '\n'
#print r1
#print '*'*10
#print t2
#print '\n'
#print r2

t = np.concatenate( (t1,t2) )
r = np.concatenate( (r1,r2) )
lamvec = np.arange(0.2,25,0.2)

for n in range(len(t)):
	print r[n]*np.cos(t[n]), r[n]*np.sin(t[n])

#print t
#print r
#print lamvec

S = []

#print t.shape
#print r.shape
#print lamvec.shape

for lam in lamvec:
	A = np.sin(2*t*k/3.)*sp.jv(2*k/3, np.sqrt(lam)*r)
	A = A.T
	[Q,R] = linalg.qr(A,mode='economic')
	#print A.shape
	#print Q.shape
	#print R.shape
	#print Q[:num_pts,:].shape
	#print '*'*5
	s = np.min(linalg.svd(Q[:num_pts-1,:])[1])
	#print s
	S.append(s)

S = np.array(S)

#for t in range(len(S)):
#	print lamvec[t], S[t]



T = range(len(S))
minima = []
for t in T[1:-1]:
	if S[t] < S[t+1] and S[t] < S[t-1]:
		minima.append(lamvec[t])

#print minima
