Send email inquiries to weswong@hsph.harvard.edu or to w2w.wong@gmail.com

Requires: ete3, scipy, numpy, collections, dill This was built using python 3.6

Multiscale Transmission model split across multiple scripts, each representing a different unit or scale

Most of the relevant scripts are in the Assets folder.

The following scripts form the hierarchy of the multiscale model: Individual.py: Individual attributes Household.py: Household attributes, tree structure, household splitting logic Node.py: Nodes used in the Household tree structure Bari.py: Bari level Village.py Village level Demographics.py: Catch all that brings the model together by defining the villages and manages all the various connections between households, baris, and vilalges

The following scripts are contain cross-file functions that are shared across the multiscale model Utility.py: collection of functions used throughout the simulation Config.py: all the parameter values are placed here

The various ipynb notebooks in the Assets folder are the exact code used for analyses and figure generation. The most relevant are: Generate_Pop_immuneProfiles.ipynb -- which takes pre-generated demographics files generated by Demographics.py and pulls out individual immunity estimates Clinical_trial.ipynb -- which takes the output from the various epidemic.py scripts and compares it to the actual clinical trial data Infection.ipynb -- contain some scripts that were used to develop the infection.py script and calibrate individual immunity Life.ipynb -- contains all the code used to calibrate the individual mortality/fertiliy statistics in Matlab Bangladesh lhs_sampler.ipynb -- contains the code needed to generate latin hyperspace cube samplings and is used as an input for epidemic_lhs_search.py

In order to run properly, synthetic populations must be created first using demographics.py

python -u Demographics.py ${index}

where index should be an integer and is only used to assign a unique identifier to the file that is generated.

The output of demographics.py, by default, is to create a population with the exact village and household size distribution as that found in Matlab (Figure 2). This requires the Taniuchi_village_start.json file to be present in the same folder. Taniuchi_village_start.json is a dictionary where: {'target_village_size': [initial bari number, # iterations]} These values were specific to the burn-in used in this study and were empirically determined through parameter sweeping the initialize_village function in Village.py

Demographics can also be created using the initialize_default_vilalges function in Demographics.py using the default parameters presented in the config.py file.

These synthetic populations form the basis of the model. To project forwards in time and simulate population after the switch, the simulate_postswitch_populations.py script should be run, after replacing line 34 with the directory path of a folder that contains populations simulated from demographics.py python simulate_postswitch_populations.py {idx} where again, the idx is an arbitrary integer used to assign a unique identifier to the file that is generated.

Once these synthetic populations are generated, the Taniuchi clinical trial can be replicated with one of three scripts: epidemic_lhs_search.py -- simulates the clinical trial from a json file specifying the parameter points to be explored. Line 29 should be replaced the file path the output of lhs_sampler.ipynb. epidemic_lhs_search.py allows for rapid exploration of parameter space using the multiscale model. To change the number of iterations run, refer to line 55 and replace it with the range that is desired

epidemic_comps_global.py -- simulates the clinical trial assuming mass action. Accepts a single parameter that indicates the value of beta_global used in the simulation. To change the number of iterations run, refer to line 51 and replace it with the range that is desired

python epidemic_lhs_search.py ${beta_global} epidemic_comps_mle.py -- simulates the clinical trial using the point estimates that best fit the clinical trial data. Can be specified to run with eithr the multiscale model or mass action model python -u epidemics_comps_mle.py ${type} ${index} where type should be "multi" for multiscale and "ma" for mass action. Note that this is case sensitive and that if the term does not match multi exactly it will run the mass action model. The index is a unique identifier that is appended to the end of the file name. To change the number of iterations run, refer to line 62 and replace it with the range that is desired

It is highly recommended that these are each run with a small number of iterations. If many simulations are needed, it is best to run them in parallel and assign each parallelized asssignment with a unique index identifier.

The outbreak.py is used to simulate mass vaccination or point importation simulations. python -u outbreak.py ${year} ${param_idx} ${param_file} ${attempt}

The first parameter, ${year} refers to the number of years post cessation and will use this value to search for an appropriate synthetic population generated by simulate_postswitch_populations.py. The directory for these files should be specified at line 46.

${param_idx} is used to determine what parameter set to use from the parameter file specified in ${param_file}. An example of param_file is in outbreak_paper_sims/s2_multi_parameters.txt This file is a json dictionary where: {'year': [[p, strain_type, year, type, iteration]]} year is the years post-Switch and the same as the first input of outbreak.py p is the proportion of individuals to be "infected" with the strain specified by strain type. Allowed strains are {'S2', and 'WPV'} When p is set to 0, a point importation outbreak is simulated when p is between 0 and 1.0, a vaccination campaign in children under 5 is performed where p represents vaccination coverage Type refers to multi or mass action (referred to as 'global' in the example file). Again, this is case sensitive and if the term multi is not used then the script will default to mass action.

iteration is one way of specifying the number of repetitions to create for this outbreak.

Due to the highly variable run times of these simulations, outbreak.py by default only runs a single iteration. To run multiple, the param_file should use the iteration parameter as a unique identifier for that particular parameter sampling. This combined with the ${attempt} parameter will allow the same code to be used to generate multiple iterations in parallel. Future iterations will likely simplify this drastically.

The various move_files.py are just different ways of consolidating simulation results once they are generated.