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

nx, ny, nz = fft_handler.get_grid_size()

# Example 1: FFT of a real-space density
rho_r = np.random.rand(nx, ny, nz) + 1j * np.random.rand(nx, ny, nz)
rho_g = fft_handler.fft(rho_r)
rho_r_back = fft_handler.ifft(rho_g)

# Example 2: FFT from atomic charge densities
rho_at_rad = [np.random.rand(100) for _ in range(10)] # 10 atoms
r = np.arange(0, 10) * 0.1 # radial grid
pos = np.random.rand(10, 3) # random positions of atoms
rho_g = fft_handler.fft_from_atom_centered(
    rho_at_rad=rho_at_rad,
    r=r,
    pos=pos
)
rho_r = fft_handler.ifft(rho_g)
```
