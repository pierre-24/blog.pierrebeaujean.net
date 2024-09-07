import numpy
import argparse
import scipy

AU_TO_KJMOL = 2625.5
AU_TO_EV = 27.212
KB = 3.166811563e-6  # Eh / K
AMU_TO_AU = scipy.constants.proton_mass / scipy.constants.electron_mass
ANG_TO_AU = 1. / scipy.constants.value('Bohr radius') * 1e-10
AU_TO_WAVENUMBER = 219474.63136320
AU_TO_KCAL = 627.5094740631


class Mode:
    """Describe one degree of freedom of the molecule.
    """

    def get_partition_function(self, T: float) -> float:
        raise NotImplementedError()

    def get_entropy(self, T: float) -> float:
        """Entropy, in Eh/K
        """
        raise NotImplementedError()

    def get_internal_energy(self, T: float) -> float:
        """Enthalpy, in Eh
        """
        raise NotImplementedError()


class QuantumHarmonicOscillator(Mode):
    """Quantum harmonic oscillator.

    Note that ZPE is included in both partition function and enthalpy.
    """
    def __init__(self, vib_energy: float):
        self.vib_energy = vib_energy

    def get_partition_function(self, T: float) -> float:
        beta = 1 / (KB * T)
        return numpy.exp(-beta * self.vib_energy / 2) / (1 - numpy.exp(-beta * self.vib_energy))

    def get_entropy(self, T: float) -> float:
        beta = 1 / (KB * T)
        S = -numpy.log(1 - numpy.exp(-beta * self.vib_energy)) + self.vib_energy * beta / (numpy.exp(beta * self.vib_energy) - 1)
        return S * KB

    def get_internal_energy(self, T: float) -> float:
        beta = 1 / (KB * T)
        E = self.vib_energy * (.5 + 1 / (numpy.exp(beta * self.vib_energy) - 1))

        return E


class HinderedRotor(Mode):
    """
    Simple hindered rotor, with potential $V(\theta) = V_0/2*(1-cos(k\theta))$, where:

    - $V_0$ is the rotational barrier,
    - $k$ is the symmetry number (*i.e.*, the number of minima/maxima).

    Some of this code is inspired by
    https://github.com/ReactionMechanismGenerator/RMG-Py/blob/main/rmgpy/statmech/torsion.pyx
    """
    def __init__(self, inertia: float, barrier: float, symmetry: int = 1, M: int = 200):
        self.inertia = inertia
        self.barrier = barrier
        self.symmetry = symmetry
        self.M = M

        self.energies = self.get_energy_levels()

    def get_hamiltonian(self) -> numpy.ndarray:
        N = 2 * self.M + 1
        H = numpy.zeros((N, N))
        col = 0

        A = 1 / (2 * self.inertia)

        for m in range(-self.M, self.M + 1):
            H[col, col] = A * m ** 2 + self.barrier / 2
            if col + self.symmetry < N:
                H[col, col + self.symmetry] -= self.barrier / 4
            if col - self.symmetry >= 0:
                H[col, col - self.symmetry] -= self.barrier / 4
            col += 1

        return H

    def get_energy_levels(self) -> numpy.ndarray:
        H = self.get_hamiltonian()
        eigenvalues, _ = numpy.linalg.eig(H)

        return numpy.sort(numpy.real(eigenvalues))

    def get_partition_function(self, T: float) -> float:
        beta = 1 / (KB * T)
        return numpy.sum(numpy.exp(-self.energies * beta)) / self.symmetry

    def get_entropy(self, T: float) -> float:
        q = self.get_partition_function(T)
        beta = 1 / (KB * T)
        exp_ = numpy.exp(-self.energies * beta)
        return KB * numpy.log(q) + numpy.sum(self.energies * exp_) / numpy.sum(exp_) / T

    def get_internal_energy(self, T: float) -> float:
        beta = 1 / (KB * T)
        exp_ = numpy.exp(-self.energies * beta)
        return numpy.sum(self.energies * exp_) / numpy.sum(exp_)


NDASHES = 22

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-I', '--inertia', type=float, help='inertia (in AMU Å²)', required=True)
    parser.add_argument('-b', '--barrier', type=float, help='barrier (in kJ mol⁻¹)', required=True)
    parser.add_argument('-f', '--frequency', type=float, help='vibrational freq (in cm⁻¹)', required=True)
    parser.add_argument('-s', '--symmetry', type=int, help='rotational symmetry', required=True)

    parser.add_argument('-T', '--temperature', type=float, default=298.15)
    args = parser.parse_args()

    # convert everything in atomic units
    E_vib_au = args.frequency / AU_TO_WAVENUMBER
    inertia_au = args.inertia * AMU_TO_AU * ANG_TO_AU ** 2
    barrier_au = args.barrier / AU_TO_KJMOL
    estimate_barrier_au = 2 * E_vib_au ** 2 * inertia_au / args.symmetry**2

    print('From frequency, estimated barrier is {:.3f} kJ mol⁻¹'.format(estimate_barrier_au * AU_TO_KJMOL))

    # correction from QHO to HR
    QHO = QuantumHarmonicOscillator(E_vib_au)
    HR = HinderedRotor(inertia_au, barrier_au, args.symmetry)

    print('Corrections at T={} K'.format(args.temperature))
    corr_enthalpy = HR.get_internal_energy(args.temperature) - QHO.get_internal_energy(args.temperature)
    corr_entropy = HR.get_entropy(args.temperature) - QHO.get_entropy(args.temperature)
    corr_mTS = - args.temperature * corr_entropy
    corr_gibbs = corr_enthalpy + corr_mTS

    print('=' * NDASHES)
    print('   U  {: 7.4f} kJ mol⁻¹'.format(corr_enthalpy * AU_TO_KJMOL))
    print('-T*S  {: 7.4f} kJ mol⁻¹    S = {:.4f} J mol⁻¹ K⁻¹'.format(
        corr_mTS * AU_TO_KJMOL,
        corr_entropy * AU_TO_KJMOL * 1e3
    ))
    print('-' * NDASHES)
    print('   A  {: 7.4f} kJ mol⁻¹'.format(corr_gibbs * AU_TO_KJMOL))
    print('=' * NDASHES)


if __name__ == '__main__':
    main()

