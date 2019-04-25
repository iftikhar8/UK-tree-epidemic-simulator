Job script reamdme:

Here, job_script.py is imported into main as a module. The simulations run as follows: each core on
the hpc (100 in total) gets a unique value of beta, each individual core then iterates through a set
of rho values [0 - 1].

From job_script.main we define the saving conventions (core_id @ param_id) where core_id is the id
for each core on the hpc and param_id is the saving label for each subsequent value of rho. An
ensemble result for beta=0.25 & rho=0.65 is saves as  0025_65_b_r.npy [b:=beta and r:=rho].

The data set is retrieved from input data based on the settings value.
Next the set of parameters is generated for each core, this depends on which type of simulation
we run, e.g. for Qro data the set of parameter values is 5.0 (being low) & 0.01 (being high) whilst
for homogeneous-random data we have [0.0, 1.0].

Finally we write all the simulation and parameter settings to the output folder as a .txt file, so
confusion over (future) analysis can be negated.
