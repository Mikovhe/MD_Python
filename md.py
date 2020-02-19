import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation


class particle:
    def __init__(self,x=0.0,y=0.0,vx=0.0,vy=0.0,force=0.0,acceleration=0.0):
        self.x = float(x)
        self.y = float(y)
        self.vx= float(vx)
        self.vy =float(vy)
        self.force = float(force)
        self.acceleration = float(acceleration)

PARTICLES = []
dt = 0.01
N = int(input("Insert number of particles\n"))
L = int(input("Insert length of box\n"))
steps = int(input("Enter number of steps\n"))
minR = 1.1
maxR = 3.2
np.random.seed(42)

def initialize():
    cvx= 0.0
    cvy=0.0
    for i in range(N):
        q = particle(float(np.random.rand(1))*L,float(np.random.rand(1))*L,float(np.random.rand(1)),float(np.random.rand(1)))
        cvx +=q.vx
        cvy +=q.vy
        PARTICLES.append(q)
        #print(center_velocity)

    for i in range(N):
        q = PARTICLES[i]
        q.vx = q.vx - cvx/N
        q.vy = q.vy-cvy/N

def verlet_algorithm():
    for i in range(N):
        p = PARTICLES[i]

        p.x = p.x + p.vx*dt + p.acceleration*0.5
        if(p.x <0.0):
            p.x  = p.x +L

        elif(p.x > L):
            p.x= p.x -L

        p.y = p.y + p.vy*dt + p.acceleration*0.5
        if(p.y<0.0):
            p.y = p.y + L

        elif(p.y>L):
            p.y = p.y - L
        PARTICLES[i] = p


def update():
    for i in range(N):
        p1 = PARTICLES[i]
        p1.force = 0.0
        for j in range(N):
            if j!=i:
                p2 = PARTICLES[j]
                r = ((p1.x-p2.x)**2 + (p1.y - p2.y)**2)**0.5
                if(r<minR):
                    f = Leonard_Force(minR)
                    p1.force += f
                    p2.force -= f
                elif(r>maxR):
                    pass
                else:
                    f=Leonard_Force(r)
                    p1.force +=f
                    p2.force -=f
        p1.acceleration = p1.force

def Leonard_Force(r):
        sigma= 1.0
        eps = 1.0
        force = -4*sigma*(-12*(eps/r)**13 + 6*(eps/r)**7)
        return force

def main():
    ims = []
    fig = plt.figure()
    initialize()
    for _ in range(steps):
        update()
        verlet_algorithm()

        X = []
        Y =[]
        for k in range(N):
            q = PARTICLES[k]
            x = q.x
            X.append([q.x])
            Y.append([q.y])
        im = plt.scatter(X,Y,color = 'black')
        ims.append([im])
    ax = fig.add_subplot(111)
    ax.set_xlim(0,L)
    ax.set_ylim(0,L)
    anim = animation.ArtistAnimation(fig,ims,interval = 800,blit = False,repeat_delay = 0)
    plt.title("Molecular dynamics simulation")
    plt.show()


main()
