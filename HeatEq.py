import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib.colors as col

def constant(a):
	return 0

def dummy(a,b):
	return 0.001*np.math.exp(-((a-1)/5)**2)

def one(a):
	return 1

class HeatEq:
	def __init__(self,a=1.0,x1=0.0,x2=1.0,t2=1.0,t_pts=100,x_pts=100,f=dummy,phi=constant,alpha=one,beta=one):
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
			self.b[0]=self.alpha(self.dt*j)
			self.b[-1]=self.beta(self.dt*j)
			for i in xrange(1,self.x_pts-1):
				self.A[i,i-1]=-coef
				self.A[i,i]=1+2*coef
				self.A[i,i+1]=-coef
				self.b[i]=self.f(j*self.dt,i*self.dx)*self.dt+(1.0-sigma)*(self.solution[j-1,i-1]+self.solution[j-1,i+1]\
					-2*self.solution[j-1,i])*self.dt*self.a*self.a/self.dx/self.dx+self.solution[j-1,i]
			self.solution[j]=np.linalg.solve(self.A,self.b)
			print self.A
			self.A.fill(0)
			self.b.fill(0)
		# self.m = np.zeros(shape=(self.x_pts, self.x_pts))
  #       self.r = np.zeros(shape=(self.x_pts,))
  #       self.solution = np.zeros(shape=(self.t_pts, self.x_pts))

  #       self.utmp = (sigma * self.dt * self.a ** 2) / (self.dx ** 2)
  #       self.dtmp = ((1 - sigma) * self.dt * self.a ** 2) / (self.dx ** 2)

  #       for i in xrange(self.x_pts):
  #           self.solution[0, i] = self.phi(i*self.dx)

  #       for t in xrange(1, self.t_pts):
  #           self.solution[t, 0] = self.alpha(self.dt*t)
  #           self.solution[t, -1] = self.beta(self.dt*t)

  #       for t in xrange(1, self.t_pts):
  #           self.m.fill(0)
  #           self.r.fill(0)
  #           self.m[0, 0] = 1
  #           self.r[0] = self.solution[t, 0]
  #           self.m[-1, -1] = 1
  #           self.r[-1] = self.solution[t, -1]
  #           self.r[0] = self.solution[t, 0]
  #           for i in xrange(1, self.x_pts - 1):
  #               self.m[i, i - 1] = -self.utmp
  #               self.m[i, i] = 1 + 2*self.utmp
  #               self.m[i, i + 1] = -self.utmp
  #               self.r[i] = self.dt*self.f(t*self.dt, i*self.dx) - self.dtmp*self.solution[t-1, i-1] + \
  #                           (1 + 2*self.dtmp)*self.solution[t-1, i] - self.dtmp*self.solution[t-1, i+1]
  #           self.solution[t] = np.linalg.solve(self.m, self.r)

	def show(self):
		def _showColor(t):
			cols=plt.pcolor((self.solution[t], self.solution[t]),\
				norm=col.Normalize(vmin=np.min(self.solution),vmax=np.max(self.solution)),cmap="plasma")
			return cols
		def _graph_animate(t):
			line.set_ydata(self.solution[t])
			return line
		fig=plt.figure(figsize=(80,5),dpi=10)
		#x = [self.x1 + i*self.dx for i in range(self.x_pts)]
		#fig, ax = plt.subplots()
		#line, = ax.plot(x, self.solution[0])
		shw=ani.FuncAnimation(fig,_showColor,frames=self.t_pts,interval=100.0/self.t_pts*100,repeat=False)
		plt.show()

	def visualize(self, type="graph"):
		def _graph_animate(t):
			line.set_ydata(self.solution[t])
			return line,
		def _pcolor_animate(t):
			cont = plt.pcolor((self.solution[t], self.solution[t]),
                              norm=col.Normalize(vmin=np.min(self.solution), vmax=np.max(self.solution)),
                              cmap="plasma")
			return cont

		x = [self.x1 + i*self.dx for i in range(self.x_pts)]
		if type == "graph":
			fig, ax = plt.subplots()
			line, = ax.plot(x, self.solution[0])
			anm = ani.FuncAnimation(fig, _graph_animate, frames=self.t_pts, interval=5000.0/self.t_pts, repeat=False)
		elif type == "pcolor":
			fig = plt.figure(figsize=(80, 5), dpi=10)
			anm = ani.FuncAnimation(fig, _pcolor_animate, frames=self.t_pts, interval=10.0/self.t_pts*100, repeat=False)
		else:
			raise Exception("Unknown plot type")

		plt.show()

a=1.0
x1=0.0
x2=1.0
t2=1000.0
t_pts=100
x_pts=100
def phi(x):
	return x/10 + np.math.exp(-(x-(x2-x1)/2)**2)
heatEq=HeatEq(a=a,x1=x1,x2=x2,t2=t2,t_pts=t_pts,x_pts=x_pts,phi=phi)
heatEq.solve(sigma=1.0/6)
heatEq.visualize(type="pcolor")