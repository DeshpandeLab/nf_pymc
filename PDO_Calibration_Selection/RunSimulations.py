from mpi4py import MPI 
import numpy as np
import sys
from HPC_exploration import model, args_run_simulations
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if __name__ == '__main__':
    PhysiCell_key, Mode, NumReplicates, Samples, Replicates = args_run_simulations(sys.argv[1:])
    NumSimulations = len(Samples)*NumReplicates
    NumSimulationsPerRank  = int(NumSimulations/size)
    mod = NumSimulations%size
    if ( mod != 0): NumSimulationsPerRank = NumSimulationsPerRank + 1 # if there is no equal split on the nodes, add a ghost simulation
    data = None

    if rank == 0:
        data = np.arange(NumSimulations) # [0,1,...,NumSimulations-1]
        if ( mod != 0):
          add = -1*np.ones((size-mod),dtype='d') # Add -1 in the ghost simulations
          data = np.concatenate((data,add),axis=None) # now len(data) % size == 0 

    recvbuf = np.empty(NumSimulationsPerRank, dtype='d')
    comm.Scatter(data, recvbuf, root=0)
    for i in range(recvbuf.shape[0]):
        if ( recvbuf[i] < 0 ): continue # If recvbuf is negative (-1) do NOT execute the model
        sampleID = int(recvbuf[i]/NumReplicates)
        if ( Mode != 'individual' ): replicateID = recvbuf[i]%NumReplicates
        else: replicateID = Replicates[sampleID]
        print('Rank: ',rank, ', Simulation: ',recvbuf[i], ', Sample: ', sampleID,', Replicate: ', replicateID)
        model(sampleID, replicateID, PhysiCell_key)
