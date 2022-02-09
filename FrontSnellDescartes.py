'''This script solves the Snell-Descarte model for seismic reflection.
It will need in input (via GUI) a .sgt file that contains the first arrival 
for different sources and receivers. 
The script will first display the first arrival in a distance-time graph where
every source (with direction +/-) is displayed independently.
Then, the user will draw the hodochrones one after the other. When all the 
hodochrones are drawn, the script will compute and show the resulting model.
For the theory behind the model: https://www.ifsttar.fr/fileadmin/user_upload/editions/lcpc/GuideTechnique/GuideTechnique-LCPC-AGAP2.pdf
'''
from matplotlib import pyplot as plt
from matplotlib import animation
from tkinter import Tk
from tkinter.filedialog import askopenfilename as askfilename
import os
import re
import numpy as np
import pygimli as pg
import time

def getH(H0:float=0.0, DipUp:float=0.0, DipDown:float=0.0, x:float=0.0):
    y0 = x*np.tan(DipUp)
    y1 = H0 + x*np.tan(DipDown)
    H = y1-y0
    return H

nbLayersDef = 2
class Snell():
    def __init__(self,nbLayers:int=nbLayersDef,V_p:np.ndarray=None,
        ThickLeft:np.ndarray=None,DipAngles:np.ndarray=None):
        '''SNELL is a class that will model the Snell-Descartes seismic refraction traveltimes
        for a given model described in input.

        INPUTS:
        -------
            - nbLayers (int) : number of layers in the model
            - V_p (np.ndarray) : array with the P-wave velocities for the different layers
            - ThickLeft (np.ndarray) : array with the thicknesses of the layer
            - DipAngles (np.ndarray): array with the dip angles
        
        All units are un SI (meters and seconds). Angles are in radians.
        '''
        self.nbLayers = nbLayers
        if V_p is not None:
            self.V_p = V_p
        else:
            self.V_p = np.linspace(1000,3000,num=self.nbLayers)
        if ThickLeft is not None:
            self.ThickLeft = ThickLeft
        else:
            self.ThickLeft = np.ones((self.nbLayers-1,))*5
        if DipAngles is not None:
            self.DipAngles = DipAngles
        else:
            self.DipAngles = np.zeros((self.nbLayers-1,))
    
    def Model(self,SourceX:float=0.0,ReceiversX:np.ndarray=np.linspace(start=0.0, stop=115.0, num=24)) -> np.ndarray:
        '''MODEL is a method that computes the forward model for the given SNELL class.

        INPUTS:
        -------
            - SourceX (float) : position of the source along X axis in m (default = 0.0)
            - RevceiversX (np.ndarray) : postions of the receivers along X axis in m
                                        (default = np.linspace(0.0,100.0,num=100))

        OUTPUT:
        -------
            - TravelTimes (np.ndarray) : travel times from sources to receivers in seconds

        '''
        # if np.count_nonzero(self.DipAngles) > 0:
        
        TravelTimes = np.zeros_like(ReceiversX)
        idx = 0
        for receiver in ReceiversX:
            # dist = abs(receiver-SourceX)
            xL = min(SourceX,receiver)
            xR = max(SourceX,receiver)
            Times = np.zeros((self.nbLayers,))
            for i in range(self.nbLayers):
                # Correction from Reynolds:
                #   - p.289 Travel-time calulations for a dipping refractor
                #   - p.285 Multilayer case
                i_crit = np.arcsin(np.divide(self.V_p[:-1],self.V_p[i])) 
                dist = abs(receiver-SourceX)
                j = 0
                DipDown = 0
                while j < i:
                    DipUp = DipDown
                    DipDown = self.DipAngles[j]
                    zL = getH(self.ThickLeft[j], DipUp, DipDown, xL)*np.cos(self.DipAngles[j])
                    zR = getH(self.ThickLeft[j], DipUp, DipDown, xR)*np.cos(self.DipAngles[j])
                    seg1 = zL/np.cos(i_crit[j])
                    seg2 = zR/np.cos(i_crit[j])
                    xL += np.sin(i_crit[j]-self.DipAngles[j])*seg1
                    xR -= np.sin(i_crit[j]-self.DipAngles[j])*seg2
                    Times[i] += (seg1+seg2)/self.V_p[j]
                    dist = np.cos(DipDown-DipUp)*dist - np.tan(i_crit[j])*(zL+zR)
                    j += 1
                if dist >= 0:
                    Times[i] += dist/self.V_p[j]
                else:
                    Times[i] = 1e20 # Absurdly high value because spacing too low to reach such deep layer.
            TravelTimes[idx] = min(Times)
            idx += 1
        return TravelTimes

    def ShowModel(self,SourceX:float=0.0,ReceiversX:np.ndarray=np.linspace(start=0.0, stop=100.0, num=101)):
        SourceX = np.asarray(SourceX)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        minX = min(SourceX.min(),ReceiversX.min())
        maxX = max(SourceX.min(),ReceiversX.max())
        meanX = (minX+maxX) / 2
        ax.plot([minX, maxX], [0, 0], 'k') # Draw the soil limit
        i = 0
        while i < self.nbLayers-1:
            minY = self.ThickLeft[:i+1].sum()
            maxY = self.ThickLeft[:i+1].sum() + np.tan(self.DipAngles[i])*(maxX-minX)
            ax.plot([minX, maxX], [minY, maxY],'k')
            ax.text(meanX, ((minY+maxY)/2*0.7), '$V_p$ = {} m/s'.format(self.V_p[i]))
            i += 1
        ax.text(meanX, max(minY,maxY)*1.1, '$V_p$ = {} m/s'.format(self.V_p[i]))
        ax.plot(ReceiversX, np.zeros_like(ReceiversX), 'kv', label='Receivers')
        ax.plot(SourceX, np.zeros_like(SourceX), 'kD',label='Source')
        ax.set_ylim([-5, max(minY,maxY)*1.5])
        ax.invert_yaxis()
        ax.set_xlabel('Distance [m]')
        ax.set_ylabel('Depth [m]')
        ax.set_title('Model display')
        ax.grid()
        fig.show()
    
    def Simulate(self, Sensors:np.ndarray, Measurements:np.ndarray, Graphs:bool=False) -> np.ndarray:
        '''SIMULATE is a method that simulates the measuremensts as performed on a field. 

        INPUTS:
        -------
            - Sensors (np.ndarray) : array containing the sensors positions
            - Measurements (np.ndarray) : array containing the configuration and measures in the form "S R T".

        OUTPUT:
        -------
            - Traveltimes (np.ndarray) : array containing the traveltimes for the configurations in Measurements
        '''
        TravelTimes = np.zeros_like(Measurements[:,-1])
        Sources = np.unique(Measurements[:,0].astype(int))
        for s in Sources:
            idx = np.where(Measurements[:,0].astype(int)==s)
            Receivers = Measurements[idx,1].astype(int)
            sX = float(Sensors[s-1,0].flatten())
            rX = Sensors[Receivers-1,0].flatten()
            t = self.Model(SourceX=sX, ReceiversX=rX)
            TravelTimes[idx] = t
        if Graphs:
            self.ShowModel(Sensors[Sources-1,0],Sensors[:,0])          
        return TravelTimes
    
    def plotHodochrones(self, Sensors:np.ndarray, Measurements:np.ndarray, Synthetic:np.ndarray=None):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        Sources = np.unique(Measurements[:,0].astype(int))
        for s in Sources:
            idx = np.where(Measurements[:,0].astype(int)==s)
            Receivers = Measurements[idx,1].astype(int)
            sX = Sensors[s-1,0].flatten()
            rX = Sensors[Receivers-1,0].flatten()
            data = Measurements[idx,2].flatten()
            ax.plot(rX,data,'d',label='Source at {} m (measured)'.format(float(sX)))
            if Synthetic is not None:
                ax.plot(rX,Synthetic[idx],':+k',label='Source at {} m (synthetic)'.format(float(sX)))
        ax.set_xlabel('Distance [m]')
        ax.set_ylabel('Travel Time [sec]')
        ax.set_title('Hodochrones')
        ax.legend(loc='upper center',
            ncol=2, fancybox=True, shadow=True)
        ax.grid()
        fig.show()

class SnellFOP(pg.core.ModellingBase): # How to build this???
    def __init__(self, sensorsX, MeasurementsArray, nlay=2, verbose=False):
        mesh = pg.meshtools.createMesh1DBlock(nlay,2) # Thicknesses, Vp and Dip
        super().__init__(mesh)
        self.x = sensorsX
        self.array = MeasurementsArray
        self.nlay = nlay
        self.nc = (nlay-1)*2 + nlay
        
    def response(self, model):
        nbLayers = self.nlay
        V_p = np.asarray(model[nbLayers-1:2*nbLayers-1])
        ThickLeft = np.asarray(model[:nbLayers-1])
        DipLeft = np.asarray(model[2*nbLayers-1:-1])
        Sensors = self.x
        Measurements = self.array
        return pg.Vector(Snell(nbLayers,V_p,ThickLeft,DipLeft).Simulate(Sensors,Measurements))
    
    def startModel(self, dataVals=None):
        ThickLeft = np.ones((self.nlay-1,))*5.0
        V_p = np.linspace(1000,3000,num=self.nlay)
        sourcesX = self.x[np.unique(self.array[:,0].astype(int))-1,0]
        if len(sourcesX) > 1:
            DipLeft = np.ones((self.nlay,))*0.001
        else:
            DipLeft = np.zeros((self.nlay,))
        modelInit = np.concatenate((ThickLeft, V_p, DipLeft))
        return pg.Vector(modelInit) # Thicknesses, Velocities, Dipping

def readSGT():
    # 1) Load a file with the first arrival:
    root = Tk()
    filename = askfilename(filetypes = (("First-Arrival", "*.sgt"), ("All types", "*.*")))
    root.destroy()
    # We retreived a first-arrival file --> geometry of the sensors + first arrivals 
    head_tail = os.path.split(filename)
    path = head_tail[0]
    filename_geom = head_tail[1]
    with open(filename) as f:
        Lines = f.read().splitlines()
    MarkerNbSensors = "# shot/geophone points"
    MarkerMeasurements = "# measurements"
    for line in Lines:
        if line.endswith(MarkerNbSensors):
            nbSensors = int(line[:-len(MarkerNbSensors)])
            Sensors = np.zeros((nbSensors,2))
            idxSensor = 0
        elif line.endswith("#x\ty"):
            pass
        elif idxSensor < nbSensors:
            Sensors[idxSensor,:] = re.split(r'\t+', line)
            idxSensor += 1
        elif line.endswith(MarkerMeasurements):
            nbMeasurements = int(line[:-len(MarkerMeasurements)])
            Measurements = np.zeros((nbMeasurements,3))
            idxMeas = 0
        elif line.endswith('#s\tg\tt'):
            pass
        elif idxMeas < nbMeasurements:
            Measurements[idxMeas,:] = re.split(r'\t+', line)
            idxMeas += 1
    return Sensors, Measurements
    
def InversionSnell(Sensors, Measurements, nlay=3):
    if np.count_nonzero(Sensors[:,1]) > 0:
        print("This dataset cannot be interpreted with the Snell refraction model!")
        print("Use PyGIMLI for the inversion of the traveltime tomography...")
        raise Exception("Dataset has topography!")
    # Get the number of sources in the system:
    Sources = Sensors[np.unique(Measurements[:,0].astype(int))-1,0]
    Receivers = Sensors[np.unique(Measurements[:,1].astype(int))-1,0]
    print("Dataset has {} sources and {} receivers.".format(len(Sources),len(Receivers)))
    ############################################################################################
    ###                                                                                      ###
    ###                                       Inversion                                      ###
    ###                                                                                      ###
    ############################################################################################
    #
    # Parameters:
    # -----------
    error = 1e-5 # AbsoluteError associated to the picking
    #
    # Inversion:
    # ----------
    data = Measurements[:,-1]
    fop = SnellFOP(Sensors, Measurements, nlay = nlay)
    # Velocities for all layers and dipping for the nlay-1 interfaces and thicknesses for all but the half-space
    regionArray = np.concatenate(((np.ones((nlay-1,))*0).astype(int), (np.ones((nlay,))*1).astype(int), (np.ones((nlay,))*2).astype(int)))
    limits = [[0.1, 100], [500, 5000], [-0.25, 0.25]]
    for nbRegion, regionType in enumerate(limits):
        fop.region(nbRegion).setLowerBound(regionType[0])
        fop.region(nbRegion).setUpperBound(regionType[1])
    inv = pg.core.Inversion(data, fop)
    inv.setLambda(0)
    inv.setRecalcJacobian(True)
    inv.setVerbose(True)
    inv.setError(error)
    inv.setDeltaPhiAbortPercent(0.001)
    model = inv.run()
    #
    # Post Processing:
    # ----------------
    V_p = np.asarray(model[nlay-1:2*nlay-1])
    ThickLeft = np.asarray(model[:nlay-1])
    DipLeft = np.asarray(model[2*nlay-1:])
    ModelContainer = Snell(nbLayers=nlay,V_p=V_p, ThickLeft=ThickLeft, DipAngles=DipLeft)
    DataSim = ModelContainer.Simulate(Sensors, Measurements)
    DataReal = Measurements[:,2]
    rms = np.sqrt(np.square(DataReal-DataSim).mean())
    print('\n\n\nResults of the inversion:\n')
    print('\t - Thicknesses [m]: {}'.format(ThickLeft))
    print('\t - Velocities [m/s]: {}'.format(V_p))
    print('\t - Dipping angles [rad]: {}\n'.format(DipLeft[:-1]))
    print('RMSE = {} seconds\n'.format(float(rms)))
    return Snell(nbLayers=nlay,V_p=V_p, ThickLeft=ThickLeft, DipAngles=DipLeft)
    
if __name__=="__main__":
    Sensors, Measurements = readSGT()
    # File has been read!
    print("First arrival times read succesfully!")
    Snell().plotHodochrones(Sensors, Measurements)
    nblay = 0
    while nblay<=1:
        try:
            nblay = int(input('\n\n\t\tInsert the number of layers for the model (nlay > 1) : \t'))
        except:
            nblay = 0
    print('\n\nBeginning inversion:')
    time.sleep(1)
    InvModel = InversionSnell(Sensors, Measurements,nblay)
    Snell(nbLayers=3,V_p=[1000,2000,4000],ThickLeft=[8,3]).Simulate(Sensors,Measurements)
    sX = Sensors[np.unique(Measurements[:,0]).astype(int)-1,0]
    rX = Sensors[:,0]
    InvModel.ShowModel(sX,rX)
    DataSimulated = InvModel.Simulate(Sensors, Measurements)
    InvModel.plotHodochrones(Sensors,Measurements,Synthetic=DataSimulated)
    plt.show()
