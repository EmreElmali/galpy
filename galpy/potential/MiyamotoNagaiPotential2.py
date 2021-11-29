import numpy
from ..util import conversion
from .Potential import Potential, kms_to_kpcGyrDecorator


class MiyamotoNagaiPotential2(Potential):

    def __init__(self, amp=1., a=1., b=0.1, normalize=False,
                 ro=None, vo=None):

        Potential.__init__(self, amp=amp, ro=ro, vo=vo, amp_units='mass')
        a = conversion.parse_length(a, ro=self._ro)
        b = conversion.parse_length(b, ro=self._ro)
        self._a = a
        self._scale = self._a
        self._b = b
        self._b2 = self._b**2.
        if normalize or \
                (isinstance(normalize, (int, float))
                 and not isinstance(normalize, bool)):
            self.normalize(normalize)
        self.hasC = False
        self.hasC_dxdv = False
        self.hasC_dens = False

    def _evaluate(self, R, z, phi=0., t=0.):
        sqrtbz = numpy.sqrt(self._b2+z**2.)
        return (-1.)/((R**2 + (self._a + sqrtbz)**2)**0.5) \
            * (1. + (self._a * (self._a + sqrtbz))/(R**2 + (self._a + sqrtbz)**2))

    def _dens(self, R, z, phi=0., t=0.):
        sqrtbz = numpy.sqrt(self._b2+z**2.)
        return (self._b2/(4*numpy.pi)) * ((3*(self._a+sqrtbz))/(sqrtbz**3*(R**2 + (self._a + sqrtbz)**2)**3.5)) \
            * (R**2*(sqrtbz**2 - self._a*sqrtbz + self._a**2) + (self._a + sqrtbz)**2 * (sqrtbz**2 + 4*self._a*sqrtbz + self._a**2))
