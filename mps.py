import numpy as np
from scipy.special import jv
from scipy import linalg

def evaluate_basis_func(angle,k,eig,point):
	r,th = point
	return jv(angle*k,np.sqrt(eig)*r)*np.sin(angle*k*np.pi)

def assemble_sample_matrix(eig, num_freq, angle, bdy_pts, int_pts):
	m = []
	for p in bdy_pts:
		row = []
		for k in range(num_freq):
			row.append(evaluate_basis_func(angle,k,eig,p))
		m.append(row)
	for q in int_pts:
		row = []
		for k in range(num_freq):
			row.append(evaluate_basis_func(angle,k,eig,p))
		m.append(row)
	return np.array(m)

def compute_sing(eig, num_freq, angle, bdy_pts, int_pts):
	sample_matrix = assemble_sample_matrix(eig, num_freq, angle,
																					bdy_pts, int_pts)
	orth = linalg.qr(sample_matrix)[0]
	sing_value = sorted(sing_vals,reverse=True)[0]
	sing_vec = sing_vecs[:,-1]
	return sing_value, sing_vec

def assemble_triangle_points(angle, r, num_pts):
	base = np.array([1,0])
	tip = np.array([r*np.cos(angle*np.pi),r*np.sin(angle*np.pi)])
	bdy_pt_list = [t*base + (1-t)*tip
	              for t in np.linspace(0,1,num=num_pts,endpoint=False)]
	int_pt_list = []
	for _ in range(num_pts):
		x = np.random.random()
		y = np.random.random()
		if x + y > 1:
			y = 1-y
			x = 1-x
		int_pt_list.append(
		                  np.array(x+base[0]*y,base[1]*y)
		                  )
	return bdy_pt_list, int_pt_list

if __name__=='__main__':
	bd,nt = assemble_triangle_points(0.5,1.0,10)
	#print bd
	#print nt
	a = assemble_sample_matrix(1.,10,0.5,bd,nt)
	#print a
	s = compute_sing(1.,10,0.5,bd,nt)[0]
	print s

	for t in np.linspace(0,15,num=1000):
		s = compute_sing(1.,10,0.5,bd,nt)[0]
		print t,s
