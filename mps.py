import numpy as np
from scipy.special import jv
from scipy import linalg
import sys
from pinp import point_in_poly
from plotting import plot_data

import matplotlib as mpl
import matplotlib.pyplot as plt

class Polygon:

	def __init__(self, bdry_pts):
		self.bdry_pts = bdry_pts

	def is_in(self,pt):
		return point_in_poly(pt[0], pt[1], self.bdry_pts)

	def bbox(self):
		"""Return xmin, xmax, ymin, ymax"""
		xlist, ylist = zip(*self.bdry_pts)
		return min(xlist), max(xlist), min(ylist), max(ylist)

	def edge_points(self, num_pts_per_edge, base=None):
		"""We assume that the bdry points are arranged in cyclic order"""
		edge_pt_list = []
		for k in range(-1,len(self.bdry_pts)-1):
			if (base is not None) and (k is not base and
			    k+1 is not base):
				strt = np.array(self.bdry_pts[k])
				fin = np.array(self.bdry_pts[k+1])
				ts = np.linspace(0., 1., num=num_pts_per_edge)
				for t in ts:
					edge_pt_list.append(strt + t*(fin - strt))
		return list(edge_pt_list)

	def find_interior_points(self, num_interior_pts):
		int_pt_list = []
		min_x, max_x, min_y, max_y = self.bbox()
		while len(int_pt_list) < num_interior_pts:
			x = (max_x-min_x)*np.random.rand() + min_x
			y = (max_y - min_y)*np.random.rand() + min_y
			if self.is_in([x,y]):
				int_pt_list.append( [x,y] )
		return list(int_pt_list)

	def translate_by(self, x, y):
		new_bdry_pts = [[p[0]+x,p[1]+y] for p in self.bdry_pts]
		return Polygon(new_bdry_pts)

	def scale_by(self, scale_factor):
		new_bdry_pts = [[scale_factor*p[0], scale_factor*p[1]]
		                 for p in self.bdry_pts]
		return Polygon(new_bdry_pts)

def to_polar(poly, vertex_index, num_pts_per_edge, num_int_pts, plot=False):
	base_pt = poly.bdry_pts[vertex_index]
	x,y = base_pt
	new_poly = poly.translate_by(-x,-y)
	edge_pts = np.array(new_poly.edge_points(num_pts_per_edge, base=vertex_index))
	int_pts = np.array(new_poly.find_interior_points(num_int_pts))
	num_int_pts = int_pts.shape
	pts = np.concatenate((edge_pts, int_pts), axis=0)
	rs = np.sqrt(pts[:,0]**2 + pts[:,1]**2)
	thetas = np.arctan2(pts[:,0], pts[:,1])
	successor = new_poly.bdry_pts[(vertex_index+1)%len(poly.bdry_pts)]
	angle_of_successor = np.arctan2(successor[1], successor[0])
	thetas = thetas - angle_of_successor
	if plot:
		ax = plt.subplot(111, projection='polar')
		ax.scatter(thetas, rs)
		ax.set_rmax(2.0)
		ax.grid(True)
		plt.show()
	return rs, thetas

def compute_angle(poly, vertex_index):
	base_pt = poly.bdry_pts[vertex_index]
	x,y = base_pt
	new_poly = poly.translate_by(-x,-y)
	successor = new_poly.bdry_pts[(vertex_index+1)%len(poly.bdry_pts)]
	predecessor = new_poly.bdry_pts[vertex_index - 1]
	angle = np.arccos( successor[0]*predecessor[0] + successor[1]*predecessor[1] )
	a = np.pi/angle
	return a

def assemble_matrix(rs, thetas, a, lam, ks, num_edges, plot=False):
	return {'matrix': evaluate_basis_func_bessel(a, lam, rs, thetas, ks)}

def u(k, a, lam, r, theta):
	return jv(a*k, np.sqrt(lam)*r)*np.sin(a*k*theta)

def old_basis_func(angle, eig, rvec, tvec, kvec):
	"""Take in index vector k, vector of rs, and vector of thetas.
	Return matrix A."""
	rows = []
	ks = kvec.flatten()
	for k in ks:
		rows.append(u(k, angle, np.sqrt(eig), rvec, tvec))
	rows = np.array(rows)
	return rows

def evaluate_basis_func_bessel(angle, eig, rvec, tvec, kvec):
	A = (jv(angle*kvec,np.sqrt(eig)*rvec)*np.sin(angle*kvec*tvec)).T
	return A

def find_sing_val(A, num_bdy_pts):
	"""Find minimal singular value of boundary part of the Q of A's QR."""
	Q,R = linalg.qr(A, mode='economic')
	#print A.shape
	#print Q.shape
	#print num_bdy_pts
	sing_val = np.min(linalg.svd(Q[:num_bdy_pts,:])[1])
	#print sing_val
	return sing_val

def compute_eigenvalues(poly):
	"""Compute Dirichlet eigenvalues of polygon"""
	return

def compute_square_eigenvalues(write_out=False):
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
	
	plt.scatter(x,y)
	plt.axes().set_aspect('equal')
	plt.show()
	
	r = np.sqrt(x**2 + y**2)
	t = np.arctan2(y,x)

	lams = np.arange(0.5, 100, 0.5)
	S = []
	Alist = []
	for lam in lams:
		A = evaluate_basis_func_bessel(a, lam, r, t, k)
		Alist.append(A)
		S.append(find_sing_val(A,4*N))#2*num_pts))

	if write_out:
		with open('shape.dat', 'w+') as f:
			for n in range(4*num_pts):
				f.write( str(r[n]*np.cos(t[n]))+' '+str(r[n]*np.sin(t[n])) + '\n')
	
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

	return lams, S, Alist

	if write_out:
		plot_data(lams, [S], ['minimum singular value'], 'minimum singular value',
	            'frequency', 'sing value', 'sing_vals.html')

def test_poly_1():
	tri = Polygon([[0.,0.],[1.,0.],[0.5,np.sqrt(3.)/2.]])
	bdy_pts = tri.edge_points(40)
	int_pts = tri.find_interior_points(250)
	pts = bdy_pts + int_pts
	with open('shape.dat','w+') as f:
		for p in pts:
			f.write(str(p[0])+' '+str(p[1])+'\n')

if __name__ == '__main__':
	compute_square_eigenvalues()
	test_poly_1()
