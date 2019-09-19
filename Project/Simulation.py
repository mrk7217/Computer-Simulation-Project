# Margaret Knoblock
# 5/4/19
# Final Project
from Body import Body
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Simulation(object):
    def __init__(self, numTimesteps, timestep, inputFileName):
        #set simulation parameters
        self.numTimesteps= numTimesteps # set by user
        self.timestep = timestep #set by user
        self.names = [] # names of bodies, will be filled through input file info
        self.allBodies = [] # will hold all body objects in system
        self.gravity = 6.67408e-11 # gravitational constant
        self.trajectories = [] # will hold the position of each body at each timestep. each body has its own row and the index in that row is the timestep
        self.patches = [] # will hold a patch object for each body
        self.numOrbits = [] # hold number of orbits for each body
        self.energyInterval = 50 #timestep interval that total energy will be calculated
        self.totalEnergies = [] # holds total energy of system at given time intervals
        self.satellite = False # boolean that is true when running simulation with a satellite
        sun_mass = 0 # variable to hold mass of sun

        #get data from input file
        file = open(inputFileName, "r") # open input file
        lines = file.readlines() # read lines of input file
        for i in range(0, len(lines), 3):
            #gets every 3rd line in input file. Format of file is Name of body \n Mass of body \n Initial position of body (if satellite then initial velocity vector values will be on this line as well)
            self.names.append(lines[i]) #name is on first line
            self.numOrbits.append(0) # initialize number of orbits for each body
            mass = float(lines[i+1]) #mass is on second line
            if i ==0: sun_mass = mass # if the sun, then save the sun mass separately
            if lines[i] == "Satellite\n": #if the object is a satellite then the third line of information will have the inital position, and initial velocity vectors
                satLine = lines[i+2].split(" ") # array to store all satellite info
                pos = float(satLine[0]) # position is first in array
            else:
                pos= float(lines[i+2]) #if not a satellite then only position is stored on this third line
            if pos == 0.0:  # if at the zero position then initial velocity is zero
                velArray = np.array([0,0])
            else: # for all other bodies
                vel = math.sqrt(self.gravity*sun_mass/pos) # calculate velocity
                if lines[i] == "Satellite\n": # in case of a satellite with its own initial velocity
                    velArray = np.array([float(satLine[1]), float(satLine[2])]) 
                    # velocity vectors for satellite are the last two values in 
                    #the array created from the third line of input for the satellite
                else:
                    velArray = np.array([0,vel]) #if not satellite make array for velocity
            tempBody = Body(mass, np.array([pos,0]), velArray, i/3) #create body for each of the inputs
            self.allBodies.append(tempBody) # add this body to the list of bodies
            
        # find initial accelerations for all bodies and store initial positions in trajectories array
        for i in range(len(self.allBodies)):
            self.allBodies[i].acc = self.calcAcc(self.allBodies[i])
            self.trajectories.append([]) #add row for each body to store that body's trajectory throughout the simulation
            self.trajectories[i].append(self.allBodies[i].pos) #add the body's initial position to the trajectories array, first element of the row just created for that body
            
    def update(self, dt):
        # updates the position, velocity, and acceleration of each body. stores new position in trajectories array
        count = 0
        for body in self.allBodies:
            newPos = body.updatePos(dt) # update position of body by one timestep
            newAcc = self.calcAcc(body) # find acceleration at time t+dt
            body.updateVel(newAcc,dt) # update velocity of body by one timestep
            self.trajectories[count].append(newPos) # add new position to trajectories array
            count += 1
    
    def updateTotalEnergy(self, timestep, f):
        # calculates total energy at a givn timestep and writes the value to a file f
        k = 0 #kinetic energy
        u = 0 #potential energy
        for i in range(len(self.allBodies)): # for each body
            speed = math.sqrt(self.allBodies[i].vel[0]**2 + self.allBodies[i].vel[1]**2)
            k += 0.5*self.allBodies[i].mass*(speed**2) # calculate kinetic energy
            for j in range(len(self.allBodies)):
                if j != i: # dont want potential energy between a body and itself
                    r_ij = math.sqrt((self.allBodies[j].pos[0] - self.allBodies[i].pos[0])**2 + (self.allBodies[j].pos[1] - self.allBodies[i].pos[1])**2)
                    u -= 0.5*(self.gravity*self.allBodies[i].mass*self.allBodies[j].mass)/r_ij # calculate potential energy
            totalEnergy = k + u # sum kinetic and potential energy
        self.totalEnergies.append(totalEnergy)
        f.write("At timestep " + str(timestep) + " total energy is: " + str(totalEnergy) + "\n")
            
        
    def calcAcc(self, i):
        # loop over list of bodies and calculate acceleration of body i due to other bodies in the system using gravitational force law. return total acceleration of body i
        totalAcc = np.array([0.0,0.0])
        for j in range(len(self.allBodies)):
            if j != i.bodyNum: # do not need to caluclate acceleration due to itself
                r_ij = self.allBodies[j].pos - i.pos # finds distance between bodies
                totalAcc += (self.allBodies[j].mass/(np.linalg.norm(r_ij)**3))*r_ij
        totalAcc *= self.gravity # multiply by gravitiational constant
        return totalAcc

    def animate(self, i):
        # loop over list of  patch objects and update their positions by one frame(timestep)
        for k in range(len(self.patches)):
            self.patches[k].center = (self.trajectories[k][i][0], self.trajectories[k][i][1])
            # resets the center of each patch to be the at the correct position for the ith timestep (k represents each planet)
        return self.patches

    def display(self):
        # this method creates the plot and planets, then runs the animation
        # set up plot
        fig = plt.figure(figsize = (8,6))
        ax = plt.axes()
        count = 0
        colors = ["y", "r", "b", "g", "c", "m"]
        # create list of patch objects and add them to plot
        for i in range(len(self.trajectories)):
            # set size of planets, sun is too big to be seen properly and satellite is too small to be seen at scale so these are special cases. the rest of planets are sized down the same
            if self.names[i] == "Sun\n":
                size = 4e10
            elif self.names[i] == "Satellite\n":
                size = self.allBodies[i].mass*200000
            else:
                size = self.allBodies[i].mass**(1.2/3)
            self.patches.append(plt.Circle((self.trajectories[i][0][0],self.trajectories[i][0][1]), size , color = colors[count], animated = True))
            ax.add_patch(self.patches[count])
            count +=1
        ax.axis('scaled')
        ax.set_xlim(-260000000000,260000000000)
        ax.set_ylim(-260000000000,260000000000)
        ax.set_xlabel("X Distance from Sun (m)")
        ax.set_ylabel("Y Distance from Sun (m)")
        ax.set_title("Inner Planet Orbital Simulation")

        # create animation
        numFrames = len(self.trajectories[0]) 
        anim = FuncAnimation(fig, self.animate, frames = numFrames, repeat = False, interval = 20, blit = True)
        plt.show() # display plot
        
    def calcOrbitPeriod(self, printResults):
        # calculates period of orbit in terms of earth years, print results if boolean parameter is true
        orbited = False
        currentTimestep = 0
        orbitPer = [] # orbit period for each planet
        orbitRatios = {} # dictionary where keys are planet names and values are the orbit periods in earth years
        for row in self.trajectories: # loop through trajectories of each planet
            while not orbited: 
                if currentTimestep >= len(row)-1:
                    # case where planet does not orbit in time given
                    orbited = True
                    orbitPer.append(0)
                elif row[currentTimestep][0]>0.0 and (row[currentTimestep+1][1]>=0.0 and row[currentTimestep][1]<0.0):
                    orbitPer.append(currentTimestep*self.timestep)
                    orbited = True
                currentTimestep+=1
            orbited= False
            currentTimestep = 0
        self.earthPeriod = orbitPer[self.names.index("Earth\n")] # holds the orbital period of the earth
        for i in range(len(self.names)):
            if self.names[i] != "Sun\n":
                #dont need orbital period in earth years for the sun
                orbitRatios[self.names[i]] = orbitPer[i]/self.earthPeriod # store key and value in dictionary
                if printResults: print(self.names[i].strip() + " has orbital period of " + str(orbitRatios[self.names[i]]) + " earth years.") #print results if printresults bool is true
    

    def run(self):
        # run solar simulation for a given number of timesteps. calculates total energy at designated intervals and writes these values to a file and array. if there is a satellite in the system then finds the time that the satellite is closest to mars and checks to see if satellite ever returns to earth
        if self.satellite:
            # variables needed to be initialized if satellite
            smallestDist = 10e500 # initial value for smallest distance from mars
            self.earthReturnTime = -1 # time for satellite to return to earth
            self.smallestTimestep = 0 # timestep where satellite is smallest distance from mars
            mars = self.trajectories[self.names.index("Mars\n")] # row for mars trajectories
            earth = self.trajectories[self.names.index("Earth\n")] # row for earth trajectories
            sat = self.trajectories[self.names.index("Satellite\n")] # row for satellite trajectories
        f = open("totalEnergy.txt", "w+") # file to write total energies to
        for i in range(self.numTimesteps): #runs simulation and stores all data in self.trajectories
            if i%self.energyInterval == 0: # calculates total energy at given interval
                self.updateTotalEnergy(i, f)
            self.update(self.timestep) # update simulation
            if self.satellite: # if running simulation with a satellite
                # calculate distance from satellite to mars
                distM = mars[i] - sat[i]
                distFromMars = math.sqrt(distM[0]**2 + distM[1]**2)
                if distFromMars < smallestDist: # if this distance is the smallest
                    smallestDist = distFromMars # set smallestDist to this distance
                    self.smallestTimestep = i #save the current timestep
                
                # calculate distance from satellite to earth
                distE = earth[i] - sat[i] 
                distFromEarth = math.sqrt(distE[0]**2 + distE[1]**2)
                if distFromEarth <= 6.37e6:
                    # if distance from center of earth is less than radius of earth then satellite will go back to earth
                    # radius of earht in meters
                    self.earthReturnTime = i # save time that this occurs
                    print("Satellite returns to Earth at timestep " + str(i))
        f.close()
        self.display() #visualizes results
        

def experiment1(sim):
    # method to run experiment 1, takes a simulation object as a parameter
    x_axis = np.arange(0, sim.numTimesteps, sim.energyInterval) # add energy in same increments as generated into totalEnergies aray
    y_axis = np.array(sim.totalEnergies) # make total energies array into np array
    plt.plot(x_axis, y_axis) #plot axes
    plt.xlabel("Timestep")
    plt.ylabel("Total Energy (J)")
    # plt.ylim((-6.22e33,-6.125e33)) #uncomment this line to see different y axes on energy graph
    plt.title("Totaly Energy vs Time in System")
    plt.tight_layout() # allows both axis titles to fit on plot
    plt.show() # show plot
    # As we can see from the plot, the difference between total energy throughout the time interval is minimal, so energy is conserved.

def experiment2(numTimesteps, timestepSize):
    # method to run experiment 2
    # need to create a different simulation object that takes in an additional body
    sim2 = Simulation(numTimesteps,timestepSize, "Experiment2-PlanetInfo.txt")
    sim2.satellite = True
    sim2.run()
    sim2.calcOrbitPeriod(False) #calculate orbital periods, dont print results
    earthDays = 365*((sim2.smallestTimestep*timestepSize)/sim2.earthPeriod) # find time to reach mars in earth days
    print("It takes " + str(earthDays) + " earth days for the satellite to reach Mars")
    if sim2.earthReturnTime != -1: # if returns to earth
        yrs = (sim2.earthReturnTime*timestepSize)/sim2.earthPeriod # calculate when in earth years this happens
        print("Satellite returns to earth in " + str(yrs) + " Earth years")
    else: # if does not return to earth
        yrs = (numTimesteps*timestepSize)/sim2.earthPeriod # calculate number of earth years the simulation runs for
        print("Satellite does not return to earth in " + str(yrs) + " Earth years")


def main():
    # creates simulation object and runs the simulation
    # runing with 1000 timesteps, each of size 100000 seconds
    sim = Simulation(1000, 100000,"PlanetInfo.txt")
    sim.run()
    sim.calcOrbitPeriod(True) # calculates the orbit period
    experiment1(sim) # runs experiment 1 on the simulation object created above
    # comment out the line above if you do not want to run experiment 1



main() # code to run basic simulation and experiment 1
# comment out the line above if you do not want to run the basic simulation and experiment 1


experiment2(1000,100000) # code to run experiment 2
# comment out the line above if you do not want to run experiment 2
