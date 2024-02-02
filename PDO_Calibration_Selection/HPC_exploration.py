import pathlib, os, sys, subprocess, time
from shutil import copyfile
import xml.etree.ElementTree as ET
import numpy as np
import random
import copy

class PhysiCell_simulations:
    def __init__(self):
        self.projName = "Template"
        self.executable='./project'
        self.configFile_ref = 'config/PhysiCell_settings.xml'
        # Parameters samples
        self.parameterSamples = np.zeros(1) # parameters samples numpy array 2D [sample_idx, parameter_idx], if np.zeros(1) is one sample no change parameters
        self.folder_configFile = None # folder to storage the config files (.xmls)
        self.configFile_name='config_S%06d_%02d.xml' # config files structure
        self.folder_outputs = "output/" # folder to storage the output folders
        self.outputs_folder = 'output_%06d_%02d/' # prefix of output folders
        self.customFile_key = None # name of parameter in xml to have a custom name
        self.customFile_value = 'File_S%06d_%02d.dat' # the structure of custom name
        # dictionary with parameters to change in the xml, parameters == None will be fill in with self.parameterSamples
        self.parameters = {}
        self.keys_variable_params = [] # list with the order of parameters to change. The parameters that will change is None.
        # HPC parameters
        self.slurm_account_name = None # account name for cluster simulation
        self.email = None # receive email if the job fail
        self.numReplicates = None # number of replicates for each simualtion
        self.omp_num_threads = None # number of threads omp for PhysiCell simulation
        self.numNodes = None # number of nodes
        self.memory = None # memory allocation to simulations
        self.numTaskPerNode = None # number of processors on each node - avaliability: max(omp_num_threads * numTaskPerNode)
        self.time = "0-00:00:00" # format of time day-hour:minute:second
    @property # decorator
    def parameters(self):
        return self._parameters
    @parameters.setter 
    def parameters(self, dicIn):
        self._parameters = dicIn
        for key, value in self.parameters.items():
            if value is None: self.keys_variable_params.append(key)
    def get_configFilePath(self,sampleID, replicateID):
        if (self.folder_configFile): 
            os.makedirs(os.path.dirname(self.folder_configFile), exist_ok=True)
            return self.folder_configFile+self.configFile_name%(sampleID,replicateID)
        else:
            folder = self.get_outputPath(sampleID, replicateID)
            return folder+self.configFile_name%(sampleID,replicateID)
    def get_outputPath(self,sampleID, replicateID):
        folder = self.folder_outputs+self.outputs_folder%(sampleID,replicateID)
        os.makedirs(os.path.dirname(folder), exist_ok=True)
        return folder
    def info(self):
        print(f"\nProject name: {self.projName}\n Executable: {self.executable}\n Config. file of reference: {self.configFile_ref}\n Parameters: {self.keys_variable_params}\n Parameters samples shape: {self.parameterSamples.shape}\n Folder to save config. files: {self.folder_configFile}\n Folder to save output folders: {self.folder_outputs}")
    def createXMLs(self): # Give a array with parameters samples generate the xml files for each simulation
        for sampleID in range(self.parameterSamples.shape[0]):
            for replicateID in range(self.numReplicates):
                ConfigFile = self.get_configFilePath(sampleID,replicateID)
                if (self.folder_outputs): self.parameters['.//save/folder'] = self.get_outputPath(sampleID, replicateID) # else save in folder of reference config file (util if there is a custom type of output)
                if (self.customFile_key): self.parameters[ self.customFile_key] = self.customFile_value%(sampleID,replicateID)
                self.parameters['.//omp_num_threads'] = self.omp_num_threads # number of threads omp for PhysiCell simulation
                self.parameters['.//random_seed'] = random.randint(0,4294967295) # random seed for each simulation
                # update the values of parameter from None to the sampled
                for idx, param_key in enumerate(self.keys_variable_params): self.parameters[param_key] = self.parameterSamples[sampleID, idx]
                generate_xml_file(pathlib.Path(self.configFile_ref), pathlib.Path(ConfigFile), self.parameters)
    def create_JOB(self, ID_Job, args): # Create a slurm script job
        print("Args: ",args_run_simulations(args)) # Checking args
        original_stdout = sys.stdout
        sys.stdout = open('job_'+str("%06d"%ID_Job)+'.sh','w')
        print ("#!/bin/bash\n")
        if (self.email):
            print ("#SBATCH --mail-user="+self.email)
            print ("#SBATCH --mail-type=FAIL")
        print ("#SBATCH --job-name="+self.projName+str("%06d"%ID_Job))
        print ("#SBATCH -p general")
        print ("#SBATCH -o "+self.projName+"_%j.txt")
        print ("#SBATCH -e "+self.projName+"_%j.err")
        print ("#SBATCH --nodes="+str(self.numNodes))
        print ("#SBATCH --ntasks-per-node="+str(self.numTaskPerNode))
        print ("#SBATCH --cpus-per-task="+str(self.omp_num_threads))
        print ("#SBATCH --time="+self.time)
        print("#SBATCH -A "+self.slurm_account_name)
        if (self.memory):
            print ("#SBATCH --mem="+self.memory)
        print ("\n module load python/3.9.8")
        print ("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
    
        run = 'python hpc/RunSimulations.py'
        for arg in args:
            run += ' '+arg
        print ("srun --cpu-bind=sockets "+run)
    
        sys.stdout = original_stdout

PhysiCell = {} # dictionary of PhysiCell structures

# Define the PhysiCell execution            
def model(Sample, Replicate, pc_key):
    ConfigFile = PhysiCell[pc_key].get_configFilePath(Sample,Replicate)
    # Write input for simulation & execute
    callingModel = [PhysiCell[pc_key].executable, ConfigFile]
    cache = subprocess.run( callingModel,universal_newlines=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if ( cache.returncode != 0):
        print(f"Error: model output error! Sample: {Sample} Replicate {Replicate}. returned: \n{str(cache.returncode)}")
        os._exit(1)

# Give a xml input create xml file output with parameters changes (verify this function for multiple cell types)
def generate_xml_file(xml_file_in, xml_file_out, parameters):
    copyfile(xml_file_in,xml_file_out)
    tree = ET.parse(xml_file_out)
    xml_root = tree.getroot()
    # key cell cycle example: xml_root.findall(".//*[@name='CD8 Tcell']/phenotype/cycle/phase_transition_rates/rate[4]")
    # key substrates example: xml_root.findall(".//*[@name='TNF']/physical_parameter_set/diffusion_coefficient")
    # key parameter example: xml_root.findall(".//random_seed")
    # Loop in all parameters
    for key in parameters.keys():
        val = parameters[key]
        # print(key,val)
        elem = xml_root.findall(key)
        if (len(elem) != 1):
            sys.exit(f"Error: multiples occurrences or none found to this key: {key}, occurences: {[pos.text for pos in elem]}")
        elem[0].text = str(val)
    tree.write(xml_file_out)

# Check if the arguments are compatible with PhysiCell simulations structure
def args_run_simulations(args):
    pc_key = args[0] 
    mode = args[1]
    if (pc_key not in PhysiCell): 
        print(f"Error: The key {pc_key} is not in the dictionary of PhysiCell structures!\n Please provide a valid key in {PhysiCell.keys()}.")
        sys.exit(1)
    try:
        if( mode == 'sequential' ): # [PhysiCell_key] sequential [Initial ID sample] [Final ID sample]
            num_replicates = PhysiCell[pc_key].numReplicates
            initial_sample = int(args[2])
            final_sample = int(args[3])
            if ( (initial_sample >= PhysiCell[pc_key].parameterSamples.shape[0]) or (final_sample > PhysiCell[pc_key].parameterSamples.shape[0]) or ( initial_sample >= final_sample) or  (initial_sample < 0) ):
                print(f"Error: Sample range unexpected! Interval [{initial_sample} , {final_sample}) doesn't match with PhysiCell['{pc_key}'] samples.")
                sys.exit(1)
            Samples = np.arange(initial_sample,final_sample)# initial_sample, initial_sample+1, ..., final_sample-1
            Replicates = np.arange(num_replicates) # 0,1,...,num_replicates-1
        elif ( mode == 'samples' ): # [PhysiCell_key] samples [ ID_sample_1  ID_sample_2 ... ID_sample_n]
            num_replicates = PhysiCell[pc_key].numReplicates
            num_samples = len(args)-2
            Samples = np.zeros(shape=(num_samples), dtype='d')
            Replicates = np.arange(num_replicates) # 0,1,...,num_replicates-1
            for sample_index in range(len(Samples)):
                sample_id = int(args[sample_index+2])
                if ( (sample_id >= PhysiCell[pc_key].parameterSamples.shape[0]) or  (sample_id < 0) ):
                    print(f"Error: Sample {sample_id} unexpected! It doesn't match with PhysiCell['{pc_key}'] samples size {PhysiCell[pc_key].parameterSamples.shape[0]}.")
                    sys.exit(1)
                Samples[sample_index] = sample_id # samples starts in args[3]
        elif ( mode == 'individual' ): # [PhysiCell_key] individual [ ID_sample_1  Replicate_ID_sample_1 ... ID_sample_n Replicate_ID_sample_n]
            num_replicates = 1
            num_samples = (len(args)-2)*0.5 # pair sample and relicate  
            if ( num_samples % 1 != 0):
                print("Error: Number of args for individual runs!\n Please provide the sample and replicate pair(s).")
                sys.exit(1)
            num_samples = int(num_samples)
            Samples = np.zeros(shape=(num_samples), dtype='d')
            Replicates = np.zeros(shape=(num_samples), dtype='d')
            for sample_index in range(len(Samples)):
                sample_id = int(args[2*sample_index+2]) # args[2], args[4], args[6] ...
                replicate_id = int(args[2*sample_index+3]) # args[3], args[5], args[7] ..
                if ( (sample_id >= PhysiCell[pc_key].parameterSamples.shape[0]) or (sample_id < 0) or (replicate_id < 0) or (replicate_id >= PhysiCell[pc_key].numReplicates) ):
                    print(f"Error: Sample: {sample_id} and replicate: {replicate_id} unexpected! It doesn't match with PhysiCell['{pc_key}'].")
                    sys.exit(1)
                Samples[sample_index] =  sample_id
                Replicates[sample_index] =  replicate_id
        else:
            print(f"Error: The mode {mode} unexpected!\n Please provide a valid mode: sequential, samples, or individual.")
            sys.exit(1)
    except:
        print("Error: Unexpected syntax!")
        print("option 1)  [PhysiCell_key] sequential [Initial ID sample] [Final ID sample]  # it doesn't include Final ID sample")
        print("option 2)  [PhysiCell_key] samples [ ID_sample_1  ID_sample_2 ... ID_sample_n] # run all replicates for each sample")
        print("option 3)  [PhysiCell_key] individual [ ID_sample_1  Replicate_ID_sample_1 ... ID_sample_n Replicate_ID_sample_n] # run individual replicate for each sample")
        sys.exit(1)
    return pc_key, mode, num_replicates, Samples, Replicates


# PhysiCell['pc_1'] = PhysiCell_simulations()
# PhysiCell['pc_1'].projName = "PDAC_chemotherapy"
# PhysiCell['pc_1'].executable='./pancreatic'
# PhysiCell['pc_1'].configFile_ref = 'config/PhysiCell_settings.xml'
# # Parameters samples
# PhysiCell['pc_1'].parameterSamples = np.array([[300,480,240,60,516,5.31667e-05]])
# PhysiCell['pc_1'].configFile_name='config_S%06d_%02d.xml'
# PhysiCell['pc_1'].folder_outputs = "output/"
# PhysiCell['pc_1'].outputs_folder = 'Sim_%06d_%02d/'
# PhysiCell['pc_1'].customFile_key = './/population_fileName'
# PhysiCell['pc_1'].customFile_value = 'PopCells_S%06d_%02d.dat'
# PhysiCell['pc_1'].parameters = {'.//full_data/enable': 'false', './/SVG/enable': 'false',".//*[@name='cancer']/phenotype/cycle/phase_durations/duration[1]": None,".//*[@name='cancer']/phenotype/cycle/phase_durations/duration[2]": None,".//*[@name='cancer']/phenotype/cycle/phase_durations/duration[3]": None,".//*[@name='cancer']/phenotype/cycle/phase_durations/duration[4]": None, ".//*[@name='cancer']/phenotype/death/*[@name='apoptosis']/phase_durations/duration": None, ".//*[@name='cancer']/phenotype/death/*[@name='apoptosis']/death_rate": None}
# # HPC parameters
# PhysiCell['pc_1'].slurm_account_name = "r00241" # account name for cluster simulation
# PhysiCell['pc_1'].email = "hlimadar@iu.edu" # receive email if the job fail
# PhysiCell['pc_1'].numReplicates = 10 # number of replicates for each simualtion
# PhysiCell['pc_1'].omp_num_threads = 8 # number of threads omp for PhysiCell simulation
# PhysiCell['pc_1'].numNodes = 100 # number of nodes
# PhysiCell['pc_1'].numTaskPerNode = 3 # number of processors on each node - avaliability: max(omp_num_threads * numTaskPerNode)
# PhysiCell['pc_1'].time = "1-00:00:00" # format of time day-hour:minute:second
#
# if __name__ == '__main__':
#     PhysiCell['pc_1'].info()
#     PhysiCell['pc_1'].createXMLs()
#     PhysiCell['pc_1'].create_JOB(0, ['pc_1',"sequential", "0", "1"])
#     PhysiCell['pc_1'].create_JOB(1, ['pc_1',"samples", "0"])
#     PhysiCell['pc_1'].create_JOB(2, ['pc_1', "individual", "0", "1", "0", "9"])
