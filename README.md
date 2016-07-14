# diffusion
Simulation of Diffusing Ligand

The ligand diffusion is modeled using the diffusion equation for the ligand concentration.

Diffusion is modeled using the diffusion equation + noise. In two dimensions (2D lattice) at each time-step the chemical concentration at lattice site `(i,j)` is updated using the following algorithm:

```
c(i,j,t+dt) = ( d*( c(i+1,j,t) + c(i-1,j,t) + c(i,j+1,t) + c(i,j-1,t) - 4*c(i,j,t)) + eta )*dt + c(i,j,t)
```
For a three dimensional system, the concentration is updated at each lattice site `(i,j,k)` with a similar algorithm. Nearest neighbors in the third dimension are added to the diffusion term and `c(i,j,k,t)` is multiplied by a factor of 6 instead of 4.

## `plr2D` directory

Model cell polarization on a 2D lattice. Cells are polarized based on the chemical concentration in their environment.

Compiling:
```
$ gfortran -c mod1* mod2*; rm *.o; gfortran plrMW2d.f90 mod*
```
