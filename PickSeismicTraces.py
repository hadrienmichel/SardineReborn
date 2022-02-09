'''This script helps to manually pick the first arrival of the seismic wave
on 
'''
from obspy import read # conda install obspy
from matplotlib import pyplot as plt
from matplotlib import animation
from tkinter import Tk
from tkinter.filedialog import askopenfilename as askfilename
from copy import deepcopy
import os
import re
import numpy as np
import time

if __name__=="__main__": # Only execute the script if called directly (it is not a function nor a module)
    root = Tk()
    filename = askfilename(filetypes = (("Geometry", "*.geometry"), ("All types", "*.*")))
    root.destroy()
    # We retreived a geometry file --> geometry of the array + position of the source(s) + names of the files (must be in the same folder) 
    head_tail = os.path.split(filename)
    path = head_tail[0]
    filename_geom = head_tail[1]
    nameSave = filename_geom[:-9]
    # Retreive the datafiles names:
    SEG2Files = []
    SourcePosition = []
    ReceiversPosition = []
    sources = False
    receivers = False
    with open(filename) as f:
        Lines = f.read().splitlines()
    for line in Lines:
        if line.startswith("SOURCES"):
            sources = True
            receivers = False
        elif line.startswith("RECEIVERS"):
            sources = False
            receivers = True
        else:
            if sources:
                tmp = re.split(r'\t+', line)
                name = tmp[0]
                CurrSource = [float(i) for i in tmp[1:]]
                SEG2Files.append(name)
                SourcePosition.append(CurrSource)
            elif receivers:
                CurrReceiver = [float(i) for i in re.split(r'\t+', line)]
                ReceiversPosition.append(CurrReceiver)
    # Check if Sources in List of Receivers --> Constitute Sensors array for output file:
    Sensors = deepcopy(ReceiversPosition)
    SourcesId = np.zeros((len(SourcePosition),))
    for source in SourcePosition:
        if not(ReceiversPosition.count(source) == 1): # If the source is not in the receivers array
            Sensors.append(source)
        SourcesId = Sensors.index(source) # We store the position of the Id position of the current source in the sensor array
    PicksSave = [] # Storing the picks and their sensor Ids
    indexFile = 0
    for name  in SEG2Files:
        _, ext = os.path.splitext(name)
        if ext == '.segy' or ext == '.sgy':
            st = read(os.path.join(path,name), 'SEGY')
            beginTime = 0
        elif ext == '.seg2' or ext == '.sg2':
            st = read(os.path.join(path,name), 'SEG2')
            beginTime = float(st[0].stats.seg2["DELAY"])
        else:
            st = read(os.path.join(path,name), 'SEG2')
            beginTime = float(st[0].stats.seg2["DELAY"])
        if len(st) != len(ReceiversPosition): # If the number of geophones does not match between the loaded array and the gemoetry
            raise Exception('The file referenced in the geometry file does not match the geometry of the array!')
        deltaT = float(st[0].stats.delta)
        # Computing the xAxes values
        nbPoints = st[0].stats.npts
        timeSEG2 = np.arange(beginTime, beginTime+nbPoints*deltaT, deltaT)
        # Picking the traces:
        Picks = [None]*len(st)
        # Create the interactive Figure:
        plt.close('all')
        figureMain = plt.figure()
        figureMain.suptitle('Seismic picker (quit figure to end/save)')
        axMain = plt.subplot2grid((2,3), (0,0), rowspan=2, colspan=2, fig=figureMain)
        axZoom = plt.subplot2grid((2,3), (0,2), rowspan=2, fig=figureMain)
        # Parameters for interactivity:
        MousePosition = 0
        currSelect = 0
        changedSelect = True # To force plot at first
        First = True
        # Functions for reactions to keyboard/mouse events:
        def changeMouse(event):
            if event.inaxes is not None:
                global MousePosition
                MousePosition = event.xdata
            return 0

        def on_key(event):
            global changedSelect, currSelect
            if event.key == "up": # Change i += 1
                currSelect += 1
                if currSelect >= len(st):
                    currSelect = 0
                changedSelect = True
            elif event.key == "down": # Change i -= 1
                currSelect -= 1
                if currSelect < 0:
                    currSelect = len(st)-1
                changedSelect = True    
            return 0

        # Excluding zoom solution adapted from: https://stackoverflow.com/questions/48446351/distinguish-button-press-event-from-drag-and-zoom-clicks-in-matplotlib
        MAX_CLICK_LENGTH = 0.5 # Max time to consider a click and not a drag.
        timeOnClick = 0

        def on_press(event):
            global timeOnClick
            timeOnClick = time.time()

        def on_release(event):
            global Picks, changedSelect, timeOnClick# Only clicks inside this axis are valid.
            if event.button == 1 and ((time.time() - timeOnClick) < MAX_CLICK_LENGTH): # If left click and not dragging accross the pannel
                Picks[currSelect] = MousePosition
                changedSelect = True
            return 0

        def AnimationZoom(i):
            global changedSelect, First
            # Change plot to go at the correct position:
            axZoom.clear()
            idx = np.greater_equal(timeSEG2,MousePosition-150*deltaT) & np.less_equal(timeSEG2,MousePosition+150*deltaT)
            timeZoom = timeSEG2[idx]
            axZoom.plot(timeZoom,st[currSelect].data[idx],color='k')
            axZoom.set_xlim(left=MousePosition-150*deltaT,right=MousePosition+150*deltaT)
            axZoom.autoscale(axis='y')
            z = axZoom.get_ylim()
            axZoom.plot([MousePosition, MousePosition],z,color='r')
            if Picks[currSelect] is not None:
                axZoom.plot([Picks[currSelect], Picks[currSelect]], z,color='g')
            axZoom.set_frame_on(False)
            axZoom.tick_params(
                axis='both',       # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                left=False,
                right=False,
                labelleft=False,
                labelbottom=False) # labels along the bottom edge are off
            if changedSelect:
                # Change red graph + pick
                if not(First):
                    limitsY = axMain.get_ylim()
                    limitsX = axMain.get_xlim()
                axMain.clear()
                i = 0
                for tr in st:
                    data = tr.data 
                    data = data/(max(data)-min(data))+i
                    if i == currSelect:
                        axMain.plot(timeSEG2,data,color='r')
                    else:
                        axMain.plot(timeSEG2,data,color='k')
                    if Picks[i] is not None:
                        axMain.plot([Picks[i], Picks[i]], [i-0.5, i+0.5],color='g')
                    i += 1
                if not(First):
                    axMain.set_xlim(left=limitsX[0],right=limitsX[1])
                    axMain.set_ylim(bottom=limitsY[0],top=limitsY[1])
                First=False
                changedSelect = False
            return 0

        plt.connect('motion_notify_event', changeMouse)
        plt.connect('key_press_event', on_key)
        plt.connect('button_press_event',on_press)
        plt.connect('button_release_event',on_release)
        ani = animation.FuncAnimation(figureMain,AnimationZoom,interval=1)

        plt.show()
        # Formatting for saving with multiple files
        for i in range(len(Picks)):
            if Picks[i] is not None:
                sId = Sensors.index(SourcePosition[indexFile])
                rId = Sensors.index(ReceiversPosition[i])
                t = Picks[i]
                PicksSave.append([sId, rId, t])
        indexFile += 1

    # Remove unused sensors from the list:
    UsedSensors = [False]*len(Sensors)
    for pick in PicksSave:
        UsedSensors[pick[0]] = True
        UsedSensors[pick[1]] = True
    OldId = range(len(Sensors))
    OldId = [i for i in range(len(Sensors)) if UsedSensors[i]]
    Sensors = [Sensors[i] for i in range(len(Sensors)) if UsedSensors[i]]
    NewId = range(len(Sensors))
    for pick in PicksSave:
        pick[0] = NewId[OldId.index(pick[0])]
        pick[1] = NewId[OldId.index(pick[1])]

    # Saving the picks in a sgt file (for interpretation in PyGimli)
    f = open(os.path.join(path,nameSave+'.sgt'),'w')# Create a new file called FirstArrival.sgt in the data directory
    nbSensors = len(Sensors)
    f.write('%d # shot/geophone points\n' % nbSensors)
    f.write('#x\ty\n')
    for i in range(nbSensors):
        f.write('%.2f\t%.2f\n' % (Sensors[i][0], Sensors[i][1]))
    nbMeas = len(PicksSave)
    f.write('%d # measurements\n' % nbMeas)
    f.write('#s\tg\tt\n')
    for i in range(nbMeas):
        f.write('%d\t%d\t%f\n' % (PicksSave[i][0]+1, PicksSave[i][1]+1, max(0,PicksSave[i][2])))
    f.close()