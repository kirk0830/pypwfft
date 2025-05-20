# PyPWFFT

A Python implementation of the Fast Fourier Transform (FFT) algorithm that is always performed in every periodic boundary condition (PBC) imposed Density Functional Theory (DFT) package. Presently this code is only designed for tutorial purposes and is not optimized for performance.

## Installation

```bash
# TBD
```

## Usage

```python
import numpy as np
from pypwfft import PWFFTHandler

real_space_lattice = np.array([
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0]
])

fft_handler = PWFFTHandler(
    lattice=real_space_lattice,
    lattice_unit='angstrom',
    ecutwfc=100,
    ecutwfc_unit='Ry'
)

# Example 1: FFT of a real-space density
nx, ny, nz = fft_handler.get_grid_size()
rho_r = np.random.rand(nx, ny, nz) + 1j * np.random.rand(nx, ny, nz)
rho_g = fft_handler.fft(rho_r)
rho_r_back = fft_handler.ifft(rho_g)

# Example 2: FFT from atomic charge densities
nat = 10 # number of atoms
rho_at_rad = [np.random.rand(100) for _ in range(nat)] # 10 atoms
r = np.arange(0, 10) * 0.1 # radial grid
pos = np.random.rand(nat, 3) # random positions of atoms
rho_g = fft_handler.fft_from_atom_centered_radials(
    at_rad=rho_at_rad,
    r=r,
    pos=pos
)
rho_r = fft_handler.ifft(rho_g)

# Example 3: FFT from atomic orbitals and density matrix (R)
nR = 10 # number of supercells
lattice = np.array([
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0]
]) # lattice vectors
iR = np.random.randint(nR, 3) # random supercell indices
nat = 10 # number of atoms
ao_rad = [np.random.rand(100) for _ in range(nat)] # 10 atoms
r = np.arange(0, 10) * 0.1 # radial grid
pos = np.random.rand(nat, 3) # random positions of atoms
denmat_r = [       np.random.rand(nat, nat) 
            + 1j * np.random.rand(nat, nat) 
            for _ in range(nR)] # density matrix
rho_g = fft_handler.fft_from_atom_centered_radials(
    data=ao_rad,
    r=r,
    pos=pos,
    iR=iR,
    denmat=denmat_r
)
rho_r = fft_handler.ifft(rho_g)

# Example 4: solve Poisson equation
nx, ny, nz = fft_handler.get_grid_size()
rho_r = np.random.rand(nx, ny, nz) + 1j * np.random.rand(nx, ny, nz)
rho_g = fft_handler.fft(rho_r)

eps = 1.0 # dielectric constant
v_g = -4*np.pi / eps * rho_g / (fft_handler.gnorm**2 + 1e-10)
v_r = fft_handler.ifft(v_g)
```
