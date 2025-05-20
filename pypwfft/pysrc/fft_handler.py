import unittest
import itertools as it
import numpy as np

def _convert_angstrom_bohr(l, reverse=False):
    '''from angstrom to bohr or vice versa'''
    return l * 1.889725989 if not reverse else l * 0.52917721092

def _convert_rydberg_ev(e, reverse=False):
    '''from Rydberg to eV or vice versa'''
    return e * 13.605698066 if not reverse else e / 13.605698066

def _reciprocal_space_vectors(Rx, Ry, Rz):
    '''returns the reciprocal space lattice vectors'''
    Rmat = np.array([Rx, Ry, Rz])
    Gmat = np.linalg.inv(Rmat).T * (2 * np.pi)
    return Gmat[0], Gmat[1], Gmat[2]

def _search_the_nearest_primitive_number(n, dividee=[2, 3, 5], larger=True):
    '''
    search the nearest primitive number of n that can be divided by
    the numbers in dividee, if larger is True, then search the
    nearest larger primitive number, otherwise search the nearest
    smaller primitive number
    '''
    if n <= 1:
        raise ValueError("Input number must be greater than 1.")
    if larger:
        while any(n % d != 0 for d in dividee):
            n += 1
    else:
        while any(n % d != 0 for d in dividee) and n > 0:
            n -= 1
    if n <= 0:
        raise ValueError("No valid primitive number found.")
    return n

def build_fft_box_by_ecut(Rx, Ry, Rz, ecut):
    '''
    build the box for performing Fast-Fourier Transform by the
    kinetic energy cutoff, this is the method of most of plane
    -wave based DFT codes, such as VASP, Quantum ESPRESSO, etc.
    With this method, you will obtained a uniformed grid in
    reciprocal space.
    
    Parameters
    ----------
    lx : np.ndarray
        lattice vector in x direction, should be a 3D vector,
        unit in Angstrom
    ly : np.ndarray
        similar to lx, but in y direction
    lz : np.ndarray
        similar to lx, but in z direction
    ecut : float
        energy cutoff of planewave for defining the fineness
        of the grid, assumed to be in Rydberg
    
    Returns
    -------
    '''
    # convert to Bohr
    ax, ay, az = map(_convert_angstrom_bohr, (Rx, Ry, Rz))
    # then tpi/bohr
    bx, by, bz = _reciprocal_space_vectors(ax, ay, az)
    
    def ekin(bx, by, bz, nbx, nby, nbz):
        '''calculate the kinetic energy of one planewave in unit
        of Rydberg'''
        return np.linalg.norm(bx * nbx + by * nby + bz * nbz) ** 2 / 2

    nbxmax, nbymax, nbzmax = map(lambda x: int(ecut // (np.linalg.norm(x) ** 2 / 2)),
                                 (bx, by, bz))
    nbxmax, nbymax, nbzmax = map(_search_the_nearest_primitive_number, 
                                 (nbxmax, nbymax, nbzmax))
    ig = [np.array([i, j, k]) for i, j, k in it.product(range(-nbxmax, nbxmax + 1),
                                                        range(-nbymax, nbymax + 1),
                                                        range(-nbzmax, nbzmax + 1))
          if ekin(bx, by, bz, i, j, k) <= ecut]
    igg = [np.linalg.norm(ig) for ig in ig]
    igg_uniq = np.unique(igg)
    # ig are the reciprocal space vectors that are within the
    # kinetic energy cutoff
    
    nx, ny, nz = map(lambda x: 2 * x + 1, (nbxmax, nbymax, nbzmax))
    nbx, nby, nbz = nx, ny, nz
    nxyz, nrxx = nx * ny * nz, nx * ny * nz
    # because we do not distribute the realspace grid, nrxx is the whole number of
    # grid points in the real space, equals to nxyz
    return {'nx': nx, 'ny': ny, 'nz': nz,
            'nbx': nbx, 'nby': nby, 'nbz': nbz,
            'nxyz': nxyz, 'nrxx': nrxx,
            'bx': bx, 'by': by, 'bz': bz,
            'ig': ig, 'igg_uniq': igg_uniq}

class TestPWFFTHandlerFreeFunctions(unittest.TestCase):
    
    def test_reciprocal_space_vectors(self):
        sc_lat = np.eye(3)
        G1x, G1y, G1z = _reciprocal_space_vectors(sc_lat[0], sc_lat[1], sc_lat[2])
        recip_sc_lat = np.array([G1x, G1y, G1z])
        self.assertTrue(np.allclose(recip_sc_lat @ sc_lat, np.eye(3) * (2 * np.pi)))
        
        fcc_lat = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]])
        G2x, G2y, G2z = _reciprocal_space_vectors(fcc_lat[0], fcc_lat[1], fcc_lat[2])
        recip_fcc_lat = np.array([G2x, G2y, G2z])
        self.assertTrue(np.allclose(recip_fcc_lat @ fcc_lat, np.eye(3) * (2 * np.pi)))

    def test_search_the_nearest_primitive_number(self):
        self.assertEqual(_search_the_nearest_primitive_number(10), 30)
        self.assertRaises(ValueError, _search_the_nearest_primitive_number, 1)
        self.assertEqual(_search_the_nearest_primitive_number(11), 30)

class PWFFTHandler:
    
    def __init__(self, 
                 lattice,
                 ecutwfc,
                 lattice_unit='angstrom',
                 ecutwfc_unit='Ry'):
        pass
    
    def get_grid_size(self):
        pass
    
    def fft(self, data):
        pass
    
    def ifft(self, data):
        pass
    
    def fft_from_atom_centered_radials(self,
                                       data,
                                       r,
                                       pos,
                                       iR=None,
                                       denmat=None):
        pass

if __name__ == "__main__":
    unittest.main()
