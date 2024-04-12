import subprocess
from time import sleep
import numpy as np

# Example script to automatically run the simulations for specified conditions
# Before running this script, make sure to compile the Rust code (run
# `cargo build --release` in the terminal)

# Parameters:
parallel_processes = 4          # number of parallel processes to run

landscapes = 100                # number of independent landscapes to run

replicates = 50                 # number of replicates to run for each landscape

mutation_rate_list = [0.001]    # list of mutation rates

size_list = ['1000']            # list of population sizes

resources = ['1', '1']          # resources provided to the population (this 
                                # should match the number of resources in the 
                                # simulation program)

a_list = [0.9]                  # epistasis parameter for the RMF model with the
                                # following limits:
                                #   - a = 1 (fully additive landscape)
                                #   - a = 0 (maximally epistatic landscape)

sigma  = 0.1                    # fitness effect sizes

correlations = [-0.9, 0., 0.9]  # correlations between the traits

models = ['--null', '--full']   # specifies if to run the null model (static 
                                # fitness landscape with no ecological 
                                # interactions) or the full model (dynamic
                                # fitness landscape with ecological interactions)

create_landscapes = False       # True: generates the fitness landscapes
                                # False: runs the simulations on a list of 
                                # existing fitness landscapes



process_list = []

if create_landscapes:
    for correlation in correlations:
        for a in a_list:
            for mutation_rate in mutation_rate_list:
                process = ['../simulations/target/release/create_landscape',
                            '--landscapes',    str(landscapes)]

                sigma_a = sigma * a
                sigma_b = sigma * np.sqrt(0.5*(1. - a**2))

                m = f'--rmf 0 {sigma_a} {sigma_a * correlation} {sigma_b} {sigma_b * correlation}'
                process.extend(m.split(' '))

                process_list.append(process)
else:
    for model in models:
        for correlation in correlations:
            for a in a_list:
                for mutation_rate in mutation_rate_list:
                    process = ['../simulations/target/release/ecoevo_landscapes',
                                '--landscapes',    '0', str(landscapes),
                                '--replicates',    str(replicates),
                                '--mutation_rate', str(mutation_rate)]

                    process.append('--size')
                    process.extend(size_list)

                    process.append('--resources')
                    process.extend(resources)

                    sigma_a = sigma * a
                    sigma_b = sigma * np.sqrt(0.5*(1. - a**2))

                    m = f'--rmf 0 {sigma_a} {sigma_a * correlation} {sigma_b} {sigma_b * correlation}'
                    process.extend(m.split(' '))

                    if len(model) > 0:
                        process.append(model)

                    process_list.append(process)


i = 0
active_processes = []
completed_processes = []

# runs while there are processes to dispatch or processes are still active
while i < len(process_list) or active_processes:
    # while there are free slots for active processes dispatch them
    while len(active_processes) < parallel_processes:
        if i < len(process_list):
            active_processes.append( subprocess.Popen(process_list[i]) )
            i += 1
        else:
            break

    # check for finished processes and remove them from the active ones
    for process in active_processes:
        if process.poll() is not None:
            completed_processes.append(process)
            active_processes.remove(process)

            print(f"{len(completed_processes)}/{len(process_list)} processes complete.", flush=True)
            break

    sleep(1.)
