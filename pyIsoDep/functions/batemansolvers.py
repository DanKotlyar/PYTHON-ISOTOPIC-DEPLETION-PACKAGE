"""batemansolvers
Three solvers are enabled here to solve the Bateman equations:

(1) ODEINT solver
-----------------
Integrate a system of ordinary differential equations with RK45 adaptive time
mesh scheme

(2) EXPM solver
---------------
Compute the matrix exponential using Pade approximation

(3) CRAM solver
---------------
CHBV computes the direct action of the matrix exponential on
a vector: y = exp(H)*x. It uses the partial fraction expansion of
the uniform rational Chebyshev approximation of type (14,14).
About 14-digit accuracy is expected if the matrix H is symmetric
negative definite. The algorithm may behave poorly otherwise.

See also PADM, EXPOKIT.

Roger B. Sidje (rbs@maths.uq.edu.au)
EXPOKIT: Software Package for Computing Matrix Exponentials.
ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
"""

import numpy as np
from pyIsoDep.functions.checkerrors import _ispositive
from scipy.linalg import solve as linsolver
from scipy.linalg import expm
from scipy.integrate import odeint

# -----------------------------------------------------------------------------
# Coefficients and poles of the partial fraction expansion
# -----------------------------------------------------------------------------

# Coefficients for IPF Cram 14
C14_ALPHA = np.array([
    +0.557503973136501826E+02 - 0.204295038779771857E+03j,
    -0.938666838877006739E+02 + 0.912874896775456363E+02j,
    +0.469965415550370835E+02 - 0.116167609985818103E+02j,
    -0.961424200626061065E+01 - 0.264195613880262669E+01j,
    +0.752722063978321642E+00 + 0.670367365566377770E+00j,
    -0.188781253158648576E-01 - 0.343696176445802414E-01j,
    +0.143086431411801849E-03 + 0.287221133228814096E-03j, ],
    dtype=np.complex128)

C14_THETA = np.array([
    -0.562314417475317895E+01 + 0.119406921611247440E+01j,
    -0.508934679728216110E+01 + 0.358882439228376881E+01j,
    -0.399337136365302569E+01 + 0.600483209099604664E+01j,
    -0.226978543095856366E+01 + 0.846173881758693369E+01j,
    +0.208756929753827868E+00 + 0.109912615662209418E+02j,
    +0.370327340957595652E+01 + 0.136563731924991884E+02j,
    +0.889777151877331107E+01 + 0.166309842834712071E+02j, ],
    dtype=np.complex128)

C14_ALPHA0 = 0.183216998528140087E-11


class CramSolver:
    """CRAM depletion solver that uses incomplete partial factorization

    A method that uses an incomplete partial factorization (IPF) for the
    Chebyshev Rational Approximation Method (CRAM), as described in:
    M. Pusa, "`Higher-Order Chebyshev Rational Approximation Method and
    Application to Burnup Equations
    <https://doi.org/10.13182/NSE15-26>`_," Nucl. Sci. Eng., 182:3, 297-318.

    Parameters
    ----------
    alpha : numpy.ndarray
        Complex residues of poles used in the factorization. Must be a
        vector with even number of items.
    theta : numpy.ndarray
        Complex poles. Must have an equal size as ``alpha``.
    alpha0 : float
        Limit of the approximation at infinity

    Attributes
    ----------
    alpha : numpy.ndarray
        Complex residues of poles :attr:`theta` in the incomplete partial
        factorization. Denoted as :math:`\tilde{\alpha}`
    theta : numpy.ndarray
        Complex poles :math:`\theta` of the rational approximation
    alpha0 : float
        Limit of the approximation at infinity

    """

    def __init__(self):
        """reset the number of partial factorization"""
        self.alpha = -C14_ALPHA
        self.theta = -C14_THETA
        self.alpha0 = C14_ALPHA0

    def solve(self, A, n0, dt):
        """Solve depletion equations using IPF CRAM

        Parameters
        ----------
        A : scipy.sparse.csr_matrix
            Sparse transmutation matrix ``A[j, i]`` desribing rates at
            which isotope ``i`` transmutes to isotope ``j``
        n0 : numpy.ndarray
            Initial compositions, typically given in number of atoms in some
            material or an atom density
        dt : float
            Time [s] of the specific interval to be solved

        Returns
        -------
        numpy.ndarray
            Final compositions after ``dt``

        """
        H = A * dt
        y = n0 * self.alpha0
        ident = np.eye(A.shape[0])
        for alpha, theta in zip(self.alpha, self.theta):
            y += np.real(linsolver(H - theta*ident, alpha*n0))
        y[y < 1E-25] = 0
        return y


class expmSolver:
    """Built-in expm solver that relies on the pade approximation"""

    def __init__(self):
        """reset values with a complete list of all the nuclides"""
        pass
    
    def solve(self, mtx, n0, dt):
        """Solve the exponential of a matrix"""
        n1 = np.dot(expm(mtx * dt), n0)
        return n1


class adaptiveOdeintSolver:
    
    def __dNdt(self, n0, t, idx):
        """function produces time rate of change for each isotope"""
        
        # Obtain the interpolated fission energy, xs, and transmutation mtx
        # -----------------------------------------------------------------
        fissE, sigf, transmutationmtx = self.dep._getInterpXS(t,\
            self.xsinterp)

        # flux is used directly
        # -----------------------------------------------------------------
        if not self.dep.flagPower:
            # calculate power for this step
            self.dep.power[idx] = (self.dep.flux[idx] * sigf * n0\
                * fissE * self.dep.volume).sum()

        # power is provided and needs to be converted to flux
        # -----------------------------------------------------------------
        else:
            self.dep.flux[idx] = self.dep.power[idx] / (
                        sigf * n0 * fissE * self.dep.volume).sum()

        # define the overall matrix to represent Bateman equations
        # -----------------------------------------------------------------
        mtxA = transmutationmtx*self.dep.flux[idx] + self.dep.decaymtx

        # solve and obtain the concentrations after a single depletion
        # -----------------------------------------------------------------
        dNdt = np.dot(mtxA, n0)
        
        return dNdt
    
    
    def __init__(self, dep, xsinterp, rtol=1E-10):
        """function initalized apdative time mesh odeint solver
        

        Parameters
        ----------
        dep : object
            depletion solver object.
        xsinterp : bool
            flag for cross section interpolation.
        rtol : float, optional
            relative convergence tolerance of isotopic concentration. The
            default is 1E-10.

        Returns
        -------
        None.

        """
        _ispositive(rtol, "relative convergence tolerance")
        self.dep = dep
        self.rtol = rtol
        self.xsinterp = xsinterp
    
    
    def solve(self, rtol=1.0e-10):
        """solve change in concentration with adaptive time mesh scheme"""
        for idx, dt in enumerate(self.dep.timesteps):
            self.dep.Nt[:, idx+1] = odeint(self.__dNdt,\
                tuple(self.dep.Nt[:, idx]), np.array([0,dt]), args=(idx,),\
                    rtol=rtol)[1,:]


class odeintSolver:
    """Solve using scipy odeint RK45 adaptive time mesh scheme"""
    
    def __init__(self, rtol=1E-10):
        """function initalizes odeint solver
        

        Parameters
        ----------
        rtol : float, optional
            Isotopic concentration convergence criteria, relative difference.
            The default is 1E-10.

        Returns
        -------
        None.

        """
        _ispositive(rtol, "relative convergence tolerance")
        self.rtol = rtol
    
    def __dNdt(self, n0, dt, flt_mtx):
        """function produces time rate of change for each isotope"""
        mtx = flt_mtx.reshape(int(len(flt_mtx)**0.5), int(len(flt_mtx)**0.5))
        return np.dot(mtx, n0)
    
    def solve(self, mtx, n0, dt):
        """solve change in concentration"""
        return odeint(self.__dNdt, tuple(n0), np.array([0, dt]),\
            args=(mtx.flatten(),), rtol=self.rtol)[1,:]
    