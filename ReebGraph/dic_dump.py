from multiR import *
import pickle
import os
import nibabel as nib
import numpy as np
import multiprocessing as mp
from tqdm import tqdm
from multiprocessing import freeze_support

def my_function(t, I, s, eps):
    T1 = I.trajectories[s]
    T2 = I.trajectories[t]
    connectEventList, disconnectEventList = findConnectDisconnectEvents(T1, T2, eps)
    return (connectEventList, disconnectEventList, t)

if __name__ == '__main__':
    freeze_support()

    # Set number of CPU cores for multiprocessing
    num_cores = mp.cpu_count()

    trkfolderO = '/Users/mansoor/Documents/Projects/DTI-Analysis/dti-scripts/ReebGraph-1/trkfolderO/' # output folder
    trkfolderI = '/Users/mansoor/Documents/Projects/DTI-Analysis/dti-scripts/ReebGraph-1/trkfolderI/' # input folder
    eps = 4 # granularity parameter

    # Get list of track files
    trknamelist = os.listdir(trkfolderI)

    for i in range(len(trknamelist)):
        trkname = trknamelist[i]
        trkpathI = os.path.join(trkfolderI, trkname)

        # Load track file
        trk = nib.streamlines.load(trkpathI)
        streamlines = trk.streamlines
        I = create_image(streamlines, eps)

        # Initialize dictionary
        dic_T = {}

        # Initialize events for each trajectory
        for s in range(len(I.trajectories)):
            T1 = I.trajectories[s]
            dic_T[s] = {}

            # Add appear event
            k1 = 0
            e1 = Event("appear", s)
            dic_T[s][k1] = [e1]

            # Add disappear event
            k2 = len(T1.points) - 1
            e2 = Event("disappear", s)
            dic_T[s][k2] = [e2]

        # Process each trajectory
        for s in range(len(I.trajectories)):
            myList = [i for i in range(s+1, len(I.trajectories))]
            inputs = tqdm(myList)

            # Parallel processing
            pool = mp.Pool(processes=num_cores)
            processed_list = pool.starmap(my_function, [(t, I, s, eps) for t in inputs])
            pool.close()

            # Process results
            for result in processed_list:
                connectEventList, disconnectEventList, t = result

                # Process connect events
                for c in connectEventList:
                    e = Event("connect", s, t, c[0], c[1])
                    dic_T[s].setdefault(c[0], []).append(e)
                    dic_T[t].setdefault(c[1], []).append(e)

                # Process disconnect events
                for c in disconnectEventList:
                    e = Event("disconnect", s, t, c[0], c[1])
                    dic_T[s].setdefault(c[0], []).append(e)
                    dic_T[t].setdefault(c[1], []).append(e)

        # Save dictionary to file
        output_path = os.path.join(trkfolderO, f"{trkname}.pickle")
        with open(output_path, 'wb') as handle:
            pickle.dump(dic_T, handle, protocol=pickle.HIGHEST_PROTOCOL)
