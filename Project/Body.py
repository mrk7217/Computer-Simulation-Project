# Margaret Knoblock
# 5/4/19
# Final Project
import math
import numpy as np

class Body(object):
    def __init__(self, m, r, v, bodyNum):
        # initialize attributes of body object
        self.mass = m
        self.pos = r # 2 element array holding positon
        self.vel = v
        self.acc = np.array([0,0]) 
        self.bodyNum = bodyNum # unique number for each body
        self.prevAcc = np.array([0,0]) # stores acceleration of body at previous timestep

    def updatePos(self, dt):
        # update position by one timestep using Beeman algorithm, dt is the timestep
        r_new = self.pos + self.vel*dt + (1/6)*(4*self.acc - self.prevAcc)*(dt**2)
        self.pos = r_new
        return r_new
    
    def updateVel(self, a, dt):
        # update velocity by one time step using Beeman algorithm, a is the new acceleration and dt is the timestep
        v_new = self.vel + (1/6)*(2*a + 5*self.acc - self.prevAcc)*dt
        self.vel = v_new
        self.prevAcc = self.acc # set previous acceleration
        self.acc = a # set new acceleration
        return v_new
        

# will  need instance variables for mass, position, velocity, acceleration
# may need additional to calculate KE. orbital period etc.
