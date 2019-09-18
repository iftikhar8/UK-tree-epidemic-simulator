DOW: 20-06-19
Author: John Holden
email: py13jh@leeds.ac.uk

Sub-grid model directory:

----------- Parameter space -----------

This directory stores code responsible for generating an epidemiological parameter-space data in terms of three variable parameters. The p-space is informed by a lattice model located in the directory "model" as a file "SSTLM.py" in short for "Stochastic-Spatio-Lattice-Model".

1. Run "run_man.sh" from terminal. This calls mkdir.py and run_main.py. The script is filled with simulation settings, here we control if the program is executed in HPC mode or Local mode. HPC mode is designed as a simple task array to be run on the HPC-arc3.leeds.ac.uk

2. mkdir.py creates the output folders to store data. Firstly, the lattice type has its own directory created. For these sub-grid simulations I will b using the 'lattice' folder. Each simulation is saved with the date and an appropriate ending tag or simluation name. Lastly, inside each individual simulatoin run there is a directory for each metric '[distance, runtime, mortality, percolation, velocity]'. Inside each directory there is a file labelled by which core on the HPC was used.

3. main.py takes in paramete arguements from run_main.sh, these parameters then define which type of simulation is run and which parameter are used over  ensemble. main.py first calls job_script.py from the job_generator module. This sets up the respect data structures ect. main.py has two modes, either it can simulate followed by animating a single realisation of the SSTLM model, storing the output as sim_anim.mp4 or it can realise an ensemble output which constitutes phase-space. The individual realisation is not setup to work on the HPC. main.py imports SSTLM.main(), this function houses the propagatio algorithm i.e. the model. 

4. SSTLM.py is a modified percolation based model with a programable dispersal kernel such that pathogens can jump to infect non-local nearest neighbours. SSTLM.py is imported into main_SSTLM.py.From SSTLM.py  main_SSTLM.py then imports the function model.main(), this function then runs the propagation algorithm. 

