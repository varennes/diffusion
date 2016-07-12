# diffusion
Simulation of Diffusing Ligand

The ligand diffusion is modeled using the diffusion equation for the ligand concentration.

## `c2D` directory

Diffusion is modeled using the diffusion equation + noise on a 2D lattice. At each time-step the chemical concentration at lattice site `(i,j)` is updated using the following algorithm:

```
c(i,j,t+dt) = d*( c(i+1,j,t) + c(i-1,j,t) + c(i,j+1,t) + c(i,j-1,t) - 4*c(i,j,t) )*dt + c(i,j,t)
```
