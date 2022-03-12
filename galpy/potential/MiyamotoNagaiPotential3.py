import numpy
from ..util import conversion
from .Potential import Potential, kms_to_kpcGyrDecorator


class MiyamotoNagaiPotential3(Potential):

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
        sqrtbz = numpy.sqrt(self._b2+z**2.)        asqrtbz= self._a+sqrtbz
        return (-1.)/(R**2.+asqrtbz**2.)**0.5 * \
                (1. + ((self._a*asqrtbz)/(R**2.+asqrtbz**2.)) - 1/3*(self._a**2.*(R**2.-2.*asqrtbz**2.))/(R**2.+asqrtbz**2.)**2.)
    def _Rforce(self, R, z, phi=0., t=0.):        sqrtbz= numpy.sqrt(self._b2+z**2.)
        asqrtbz= self._a+sqrtbz        return -R/(R**2.+asqrtbz**2.)**(1.5) * \                (1. + ((self._a*asqrtbz)/(R**2.+asqrtbz**2.)) - 1/3*(self._a**2.*(R**2.-2.*asqrtbz**2.))/(R**2.+asqrtbz**2.)**2.) \               -1./((R**2.+asqrtbz**2.))**0.5 * \                ((-2.*R*self._a*asqrtbz)/(R**2.+asqrtbz**2.)**2. + 2.*self._a**2.*R/3. * (R**2-5.*asqrtbz**2)/(R**2.+asqrtbz**2.)**3.  )    def _zforce(self, R, z, phi=0., t=0.):        sqrtbz= numpy.sqrt(self._b2+z**2.)
        asqrtbz= self._a+sqrtbz        return (-z*asqrtbz/sqrtbz/(R**2.+asqrtbz**2.)**(1.5)) * \                (1. + ((self._a*asqrtbz)/(R**2.+asqrtbz**2.)) - 1/3*(self._a**2.*(R**2.-2.*asqrtbz**2.))/(R**2.+asqrtbz**2.)**2.) \
               -1./((R**2.+asqrtbz**2.))**0.5 * \                (self._a*z/sqrtbz * (R**2-asqrtbz**2)/(R**2.+asqrtbz**2.)**2. + (4.*self._a**2*sqrtbz*z/(3.*sqrtbz) * (2*R**2-asqrtbz**2)/(R**2.+asqrtbz**2.)**3.))
    def _dens(self, R, z, phi=0., t=0.):
        sqrtbz = numpy.sqrt(self._b2+z**2.)
        return (self._b2/(4.*numpy.pi)) * (1.)/(sqrtbz**3. * (R**2. + (self._a + sqrtbz)**2.)**4.5) \
            * (3.*R**4. * sqrtbz**3. + R**2.*(self._a + sqrtbz)**2. * (6.*sqrtbz**3. + 15.*self._a*sqrtbz**2. - 10.*self._a**2.*sqrtbz + 5.*self._a**3.) \
                + (self._a + sqrtbz)**4. * (3.*sqrtbz**3. + 15.*self._a*sqrtbz**2. + 25.*self._a**2.*sqrtbz + 5.*self._a**3.))
    def _R2deriv(self,R,z,phi=0.,t=0.):        sqrtbz= numpy.sqrt(self._b2+z**2.)
        asqrtbz= self._a+sqrtbz        return (1./(R**2.+asqrtbz**2.)**1.5 - 3.*R**2./(R**2.+asqrtbz**2.)**2.5) * \                (1. + ((self._a*asqrtbz)/(R**2.+asqrtbz**2.)) - 1/3*(self._a**2.*(R**2.-2.*asqrtbz**2.))/(R**2.+asqrtbz**2.)**2.) \               -R/(R**2.+asqrtbz**2.)**(1.5) * \                ((-4.*R*self._a*asqrtbz)/(R**2.+asqrtbz**2.)**2. + 4.*self._a**2*R/3. * (R**2.-5.*asqrtbz**2.)/(R**2.+asqrtbz**2.)**3.) \               -1./((R**2.+asqrtbz**2.))**0.5 * \                ((24.*R**2*self._a*asqrtbz - 6.*self._a*asqrtbz*(R**2.+asqrtbz**2.) + 2.*self._a**2.*(R**2.-5.*asqrtbz**2.))/(3.*(R**2.+asqrtbz**2.)**3.) + (4.*R**2*self._a**2.*(R**2.+asqrtbz**2.) - 12.*R**2.*self._a**2.*(R**2.-5.*asqrtbz**2.))/(3.*(R**2.+asqrtbz**2.)**4.))        # First term starts with the second derivative of Phi_MN1 in z; here, Bovy's implementation (as in the original "MiyamotoNagaiPotential") is being used.    # Alternatively, Michael's equivalent 2nd z-derivative can be implemented.    def _z2deriv(self,R,z,phi=0.,t=0.):        sqrtbz= numpy.sqrt(self._b2+z**2.)
        asqrtbz= self._a+sqrtbz        return ((self._a**3.*self._b2 + 
                     self._a**2.*(3.*self._b2 - 2.* z**2.)
                     *numpy.sqrt(self._b2 + z**2.)
                     + (self._b2 + R**2. - 2.*z**2.)*(self._b2 + z**2.)**1.5
                     +self._a* (3.*self._b2**2. - 4.*z**4. + self._b2*(R**2. - z**2.)))/
                     ((self._b2 + z**2.)**1.5* (R**2. + asqrtbz**2.)**2.5)) * \                        (1. + (self._a*asqrtbz)/(R**2.+asqrtbz**2.) - 1/3*self._a**2.*(R**2.-2.*asqrtbz**2.)/(R**2.+asqrtbz**2.)**2.) \                     -z*asqrtbz/sqrtbz/(R**2.+asqrtbz**2.)**1.5 * \                        (2.*self._a*z/sqrtbz * (R**2.-asqrtbz**2)/(R**2.+asqrtbz**2.)**2. +  8.*self._a**2.*asqrtbz*z/(3.*sqrtbz) * (2.*R**2.-asqrtbz**2)/(R**2.+asqrtbz**2.)**3.) \                     -1./((R**2.+asqrtbz**2.))**0.5 * \                        ((self._a*(sqrtbz**2.-z**2)/sqrtbz**3.) * (R**2.-asqrtbz**2)/(R**2.+asqrtbz**2.)**2. - (2.*self._a*asqrtbz*z**2./sqrtbz**2.) * (3.*R**2.-asqrtbz**2)/(R**2.+asqrtbz**2.)**3. + (4.*self._a**2.*sqrtbz**2.*asqrtbz + 4.*self._a**2.*sqrtbz*z**2.-4.*self._a**2.*asqrtbz*z**2.)/(3.*sqrtbz**3.) * (2.*R**2.-asqrtbz**2.)/(R**2.+asqrtbz**2.)**3. - (8.*self._a**2.*asqrtbz**2.*z**2.)/(3.*sqrtbz**2.) * (7.*R**2.-2*asqrtbz**2.)/(R**2.+asqrtbz**2.)**4.)