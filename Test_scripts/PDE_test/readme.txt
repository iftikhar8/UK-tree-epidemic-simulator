1. main(L,beta) is called from pde_model_main.py in parent directory.
2. based on the file name, phase_3d name, we generate a 3d velocity phase-space
3. each region inside the uk is mapped to a velocity-value from the phase space
	- using a cg- boxing method: define and upper and lower boundary of density
	  based on the phase space dimension.
	- if the density value inside the domain falls within this boundary, we 
	 assign it the corresponding velocity value.
4. we convert to a diffusion constant for our RD equation using : v = 2 * sqrt(u) --> u = v^2 / 4
5. the PDE model is iterated over using these values
