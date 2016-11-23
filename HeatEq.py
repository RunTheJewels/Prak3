import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib.colors as col


def constant(a):
	return 0

def dummy(a,b):
	return 0.01*np.math.exp(-((a-1)/5)**2)*(np.sin(a/np.pi)+1)

def one(a):
	return 1

class HeatEq:
	def __init__(self,a=1.0,x1=0.0,x2=1.0,t2=1.0,t_pts=100,x_pts=100,f=dummy,phi=constant,alpha=one,beta=one,type_alpha=1,type_beta=1):
		self.a=a
		self.phi=phi
		self.x1=x1
		self.x2=x2
		self.t2=t2
		self.f=f
		self.t_pts=t_pts
		self.x_pts=x_pts
		self.dx=(x2-x1)/x_pts
		self.dt=t2/t_pts
		self.alpha=alpha
		self.beta=beta
		if ((type_beta==1) or (type_beta==2)):
			self.type_beta=type_beta
		else:
			raise Exception("Type of condition unknown")
		if ((type_alpha==1) or (type_alpha==2)):
			self.type_alpha=type_alpha
		else:
			raise Exception("Type of condition unknown")

	def solve(self,sigma):
		self.solution=np.zeros((self.t_pts,self.x_pts))
		self.A=np.zeros((self.x_pts,self.x_pts))
		self.b=np.zeros((self.x_pts))
		coef=self.a*self.a/self.dx/self.dx*sigma*self.dt
		for i in xrange(self.x_pts):
			self.solution[0,i]=self.phi(i*self.dx)
		for j in xrange(1,self.t_pts):
			self.solution[j,0]=self.alpha(self.dt*j)
			self.solution[j,-1]=self.beta(self.dt*j)
		for j in xrange(1,self.t_pts):
			self.A[0,0]=1
			self.A[-1,-1]=1
			if (self.type_alpha==1):
				self.b[0]=self.alpha(self.dt*j)
			else:
				self.A[0,1]=-1
				self.b[0]=self.alpha(self.dt*j)*self.dx
			if (self.type_beta==1):
				self.b[-1]=self.beta(self.dt*j)
			else:
				self.A[-1,-2]=-1
				self.b[-1]=self.beta(self.dt*j)*self.dx
			for i in xrange(1,self.x_pts-1):
				self.A[i,i-1]=-coef
				self.A[i,i]=1+2*coef
				self.A[i,i+1]=-coef
				self.b[i]=self.f(j*self.dt,i*self.dx)*self.dt+(1.0-sigma)*(self.solution[j-1,i-1]+self.solution[j-1,i+1]\
					-2*self.solution[j-1,i])*self.dt*self.a*self.a/self.dx/self.dx+self.solution[j-1,i]
			self.solution[j]=np.linalg.solve(self.A,self.b)
			self.A.fill(0)
			self.b.fill(0)

	def show(self):
		def _graph_animate(t):
			line.set_ydata(self.solution[t])
			return line
		xs = [self.x1 + i*self.dx for i in range(self.x_pts)]
		fig, ax = plt.subplots()
		line, = ax.plot(xs, self.solution[0])
		anm = ani.FuncAnimation(fig, _graph_animate, frames=self.t_pts, interval=50000.0/self.t_pts, repeat=False)
		plt.show()


a=1.0
x1=0.0
x2=10.0
t2=80.0
t_pts=1000
x_pts=100
t_a=2
t_b=1

def phi(x):
	return np.math.exp(-(x-(x2-x1)/2)**2)

heatEq=HeatEq(a=a,x1=x1,x2=x2,t2=t2,t_pts=t_pts,x_pts=x_pts,phi=constant,type_alpha=t_a,type_beta=t_b,alpha=np.sin,beta=one)
heatEq.solve(sigma=1.0/2)
heatEq.show()