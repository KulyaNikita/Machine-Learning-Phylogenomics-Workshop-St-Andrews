import subprocess
import os
import numpy as np
import random
import sys

tmpdir   = os.getcwd() + "/"

runID = str(os.getpid()*random.randint(1,100000))

num_of_replicates = sys.argv[1]
substitution_model = sys.argv[2]
polymosim = sys.argv[3]

#global variables:
global_replicates = int(num_of_replicates)
global_replicates_short_internalbranch = int(num_of_replicates)
total_replicates = 3*global_replicates  + 3*global_replicates_short_internalbranch
global_seqlen = 100000
global_min_branchlength = 0.1
global_max_branchlength = 0.5
global_min_short_internal_branchlength = 0.001
global_max_short_internal_branchlength = 0.02
# All models
global_min_pinv   = 0
global_max_pinv   = 0.50
global_min_shape  = 0.01
global_max_shape  = 4
# K2P, F84, HKY
global_min_tstv = 1.0
global_max_tstv = 3.0
# F81, F84, HKY, GTR
global_min_basefreq = 0.2
global_max_basefreq = 0.3
# GTR
global_min_rrates_GTR = 0.1
global_max_rrates_GTR = 1.0

def find_frequency(freq_file):
    with open(freq_file,"r") as frequency_file:
        frequency_list = []
        for line in frequency_file:
            column = line.split("\t")
            if len(column[0]) == 4:
                frequency_list.append(np.float32(column[1]))
            else :
                f_e = open('file_error.txt', 'a')
                f_e.write(line)
                f_e.write('\n')
                f_e.close()
    return frequency_list


def make_tree_file_model(file_name,topology,B1,B2,B3,B4,B5,B6, seqlen, modelname):
    with open(file_name,"w") as tree_file:
        if topology == 0:
            string = "1.0  {0} {1} ((A:{2},B:{3}):{4},(C:{5},D:{6}):{7})\n".format(seqlen, modelname, B1,B2,B3,B4,B5,B6)
        if topology ==1 :
            string = "1.0  {0} {1} ((A:{2},C:{3}):{4},(B:{5},D:{6}):{7})\n".format(seqlen, modelname, B1,B2,B3,B4,B5,B6)
        if topology ==2:
            string = "1.0  {0} {1} ((A:{2},D:{3}):{4},(C:{5},B:{6}):{7})\n".format(seqlen, modelname, B1,B2,B3,B4,B5,B6)
        print(string)
        tree_file.write(string)


def make_model_file(file_name, model, alpha, pinv, tstv,basefreq_1,basefreq_2,GTR_A_C,GTR_A_G,GTR_A_T,GTR_C_G,GTR_C_T,GTR_G_T, model_name):
    with open(file_name, "w") as model_file:
        if model == 'JC':
            model_file.write("begin model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n'.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            model_file.write("end model\n")
        elif model == 'K2P':
            model_file.write("begin model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n'.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            string = 'tstv: {0}\n'.format(tstv)
            model_file.write(string)
            model_file.write("end model\n")
        elif model == 'F81':
            model_file.write("begin model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n'.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            string = 'basefreq: {0}  {1}  {2}  {3}\n'.format(basefreq_1,basefreq_2,basefreq_1,basefreq_2)
            model_file.write(string)
            model_file.write("end model\n")
        elif model == 'HKY':
            model_file.write("begin model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n'.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            string = 'tstv: {0}\n'.format(tstv)
            model_file.write(string)
            string = 'basefreq: {0}  {1}  {2}  {3}\n'.format(basefreq_1, basefreq_2, basefreq_1, basefreq_2)
            model_file.write(string)
            model_file.write("end model\n")
        elif model == 'F84':
            model_file.write("begin model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n'.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            string = 'tstv: {0}\n'.format(tstv)
            model_file.write(string)
            string = 'basefreq: {0}  {1}  {2}  {3}\n'.format(basefreq_1, basefreq_2, basefreq_1, basefreq_2)
            model_file.write(string)
            model_file.write("end model\n")
        elif model == 'GTR':
            model_file.write("begin model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n'.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            string = 'basefreq: {0}  {1}  {2}  {3}\n'.format(basefreq_1, basefreq_2, basefreq_1, basefreq_2)
            model_file.write(string)
            string = 'rrates: {} {} {} {} {} {}\n'.format(GTR_A_C,GTR_A_G,GTR_A_T,GTR_C_G,GTR_C_T,GTR_G_T)
            model_file.write(string)
            model_file.write("end model\n")
        elif model in base_models_choice_aa :
            model_file.write("begin aa-model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n '.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            model_file.write("end model\n")
            print("Model file written:")



def simulate_random_lengths_random_free_model(topology, seed, min_brlen, max_brlen, min_internal_brlen,
                                              max_internal_brlen, min_pinv, max_pinv,
                                              min_shape, max_shape, seqlen, global_min_tstv,global_max_tstv,
                                              global_min_rrates_GTR,global_max_rrates_GTR):
    model_filename = tmpdir + "Modelfile_free" + runID + "_" + str(seed) + ".txt"
    seq_filename = tmpdir + "sim-sequences" + runID + "_" + str(seed) + ".txt"
    sim_filename = tmpdir + "sim_tree" + runID + "_" + str(seed) + ".txt"


    model = substitution_model
    modelname = model + '_model'


    shape = random.uniform(min_shape, max_shape)
    pinv = random.uniform(min_pinv, max_pinv)

    r1 = random.uniform(min_brlen, max_brlen)
    r2 = random.uniform(min_brlen, max_brlen)
    r3 = random.uniform(min_internal_brlen, max_internal_brlen)
    r4 = random.uniform(min_brlen, max_brlen)
    r5 = random.uniform(min_brlen, max_brlen)
    r3 = r3 / 2
    r6 = r3
    tstv = random.uniform(global_min_tstv, global_max_tstv)
    basefreq_1 = np.random.normal(0.5,0.03,1)/2
    while basefreq_1 > 0.3 or basefreq_1 < 0.2:
        basefreq_1 = np.random.normal(0.5,0.03,1)/2
    basefreq_2 = (1 - 2*basefreq_1)/2
    basefreq_1 = float(basefreq_1)
    basefreq_2 = float(basefreq_2)
    GTR_A_C = random.uniform(global_min_rrates_GTR,global_max_rrates_GTR)
    GTR_A_G = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)
    GTR_A_T = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)
    GTR_C_G = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)
    GTR_C_T = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)
    GTR_G_T = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)

    print("Calling make_model_file: ", model_filename, " ", model, " ", shape, " ", pinv, " ", modelname)
    make_model_file(model_filename, model, shape, pinv, tstv,basefreq_1,basefreq_2,GTR_A_C,GTR_A_G,GTR_A_T,GTR_C_G,GTR_C_T,GTR_G_T, modelname)
    make_tree_file_model(sim_filename, topology, r1, r2, r3, r4, r5, r6,  seqlen, modelname)

    all_sim_filename = "all_sim_tree" + runID + ".txt"
    os.system("echo " + str(seed) + " >> " + all_sim_filename)
    os.system("cat " + sim_filename + " >> " + all_sim_filename)
    command = polymosim + " -s " + str(
        seed) + " -m " + model_filename + " -t " + sim_filename + " -f site_pattern_freq_relative_fill -n 1 1> " + seq_filename
    subprocess.run([command], shell=True)
    individual_frequencies = find_frequency(seq_filename)
    os.remove(seq_filename)
    os.remove(sim_filename)
    os.remove(model_filename)
    return individual_frequencies


#my main program starts here:
topology_list = []
topology_array = np.empty((total_replicates,1))
frequency_array = np.empty((total_replicates,256))

frequency_simulations_filename = "MyFirstSitePatternFrequencies_" + substitution_model + ".npy"
topology_simulations_filename  = "MyFirstTreeTopologies_" + substitution_model + ".npy"

for i in range(0, global_replicates):
    temporary_array = simulate_random_lengths_random_free_model(0, i,
                                                                global_min_branchlength, global_max_branchlength,
                                                                global_min_branchlength, global_max_branchlength,
                                                                global_min_pinv, global_max_pinv, global_min_shape,
                                                                global_max_shape, global_seqlen,
                                                                 global_min_tstv, global_max_tstv,
                                                                global_min_rrates_GTR, global_max_rrates_GTR)

    frequency_array[i] = temporary_array
    topologies = topology_list.append(0)
for i in range(global_replicates, 2 * global_replicates):
    temporary_array = simulate_random_lengths_random_free_model(1, i,
                                                                global_min_branchlength, global_max_branchlength,
                                                                global_min_branchlength, global_max_branchlength,
                                                                global_min_pinv, global_max_pinv, global_min_shape,
                                                                global_max_shape, global_seqlen,

                                                                global_min_tstv, global_max_tstv,
                                                                global_min_rrates_GTR, global_max_rrates_GTR)

    frequency_array[i] = temporary_array
    topologies = topology_list.append(1)
for i in range(2 * global_replicates, 3 * global_replicates):
    temporary_array = simulate_random_lengths_random_free_model(2, i,
                                                                global_min_branchlength, global_max_branchlength,
                                                                global_min_branchlength, global_max_branchlength,
                                                                global_min_pinv, global_max_pinv, global_min_shape,
                                                                global_max_shape, global_seqlen,

                                                                global_min_tstv, global_max_tstv,
                                                                global_min_rrates_GTR, global_max_rrates_GTR)
    frequency_array[i] = temporary_array
    topologies = topology_list.append(2)

for i in range(3 * global_replicates, 3 * global_replicates + global_replicates_short_internalbranch):
    temporary_array = simulate_random_lengths_random_free_model(0, i,
                                                                global_min_branchlength, global_max_branchlength,
                                                                global_min_short_internal_branchlength,
                                                                global_max_short_internal_branchlength,
                                                                global_min_pinv, global_max_pinv, global_min_shape,
                                                                global_max_shape, global_seqlen,

                                                                global_min_tstv, global_max_tstv,
                                                                global_min_rrates_GTR, global_max_rrates_GTR)

    frequency_array[i] = temporary_array
    topologies = topology_list.append(0)
for i in range(3 * global_replicates + global_replicates_short_internalbranch,
               3 * global_replicates + 2 * global_replicates_short_internalbranch):
    temporary_array = simulate_random_lengths_random_free_model(1, i,
                                                                global_min_branchlength, global_max_branchlength,
                                                                global_min_short_internal_branchlength,
                                                                global_max_short_internal_branchlength,
                                                                global_min_pinv, global_max_pinv, global_min_shape,
                                                                global_max_shape, global_seqlen,
                                                                 global_min_tstv, global_max_tstv,
                                                                global_min_rrates_GTR, global_max_rrates_GTR)

    frequency_array[i] = temporary_array
    topologies = topology_list.append(1)
for i in range(3 * global_replicates + 2 * global_replicates_short_internalbranch,
               3 * global_replicates + 3 * global_replicates_short_internalbranch):
    temporary_array = simulate_random_lengths_random_free_model(2, i,
                                                                global_min_branchlength, global_max_branchlength,
                                                                global_min_short_internal_branchlength,
                                                                global_max_short_internal_branchlength,
                                                                global_min_pinv, global_max_pinv, global_min_shape,
                                                                global_max_shape, global_seqlen,
                                                               global_min_tstv, global_max_tstv,
                                                                global_min_rrates_GTR, global_max_rrates_GTR)

    frequency_array[i] = temporary_array
    topologies = topology_list.append(2)

topology_array = topology_list

np.save(frequency_simulations_filename, frequency_array)
np.save(topology_simulations_filename, topology_array)
