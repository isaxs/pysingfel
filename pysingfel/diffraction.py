import numba as nb
from numba import jit
from scipy.interpolate import CubicSpline

from pysingfel.detector import *
from pysingfel.particle import *


def calculate_Thomson(ang):
    # Should fix this to accept angles mu and theta
    re = 2.81793870e-15  # classical electron radius (m)
    P = (1 + np.cos(ang)) / 2.
    return re**2 * P  # Thomson scattering (m^2)


def calculate_atomicFactor(particle, detector):
    """
    particle.ffTable contains atomic form factors against particle.qSample.
    Here the atomic form factors related to the wavevector of each pixel of 
    the detector is calculated through interpolation to the ffTable data.
    
    Return:
        f_hkl: atomic form factors related to each pixel of detector for each atom type.
    """
    f_hkl = np.zeros((detector.py, detector.px, particle.numAtomTypes))
    q_mod_Bragg = detector.q_mod * 1e-10/2.
    for atm in range(particle.numAtomTypes):
        cs = CubicSpline(particle.qSample, particle.ffTable[atm, :])  # Use cubic spline
        f_hkl[:, :, atm] = cs(q_mod_Bragg)  # interpolate
    return f_hkl


@jit
def Phase(atomPos,  q_xyz):
    phase = 2 * np.pi * (atomPos[0] * q_xyz[:, :, 0] + atomPos[1] * q_xyz[:, :, 1] + atomPos[2] * q_xyz[:, :, 2])
    return np.exp(1j * phase)


@jit
def cal(f_hkl, atomPos, q_xyz, xyzInd):
    F = np.zeros_like(q_xyz[:, :, 0], dtype=nb.c16)
    for atm in range(atomPos.shape[0]):
        F += Phase(atomPos[atm, :], q_xyz) * f_hkl[:, :, xyzInd[atm]]
    return np.abs(F)**2


def calculate_molecularFormFactorSq(particle, detector):
    """
    Calculate molecular form factor for each pixel of detector to get diffraction pattern. 
    Sum over all atoms with the right phase factor.
    See https://www.nature.com/article-assets/npg/srep/2016/160425/srep24791/extref/srep24791-s1.pdf
    for more details about the derivation (equ. 13 is used for calculation here).
    """
    f_hkl = calculate_atomicFactor(particle, detector)
    s = particle.SplitIdx
    xyzInd = np.zeros(s[-1], dtype=int)
    for i in range(len(s)-1):
        xyzInd[s[i]:s[i+1]] = i
    return cal(f_hkl, particle.atomPos, detector.q_xyz, xyzInd)


def calculate_compton(particle, detector):
    """
    Calculate the contribution to the diffraction pattern from compton scattering.
    """
    half_q = detector.q_mod * 1e-10/2.
    cs = CubicSpline(particle.comptonQSample, particle.sBound)
    S_bound = cs(half_q)
    Compton = S_bound + particle.nFree
    return Compton
