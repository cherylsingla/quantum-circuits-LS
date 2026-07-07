import numpy as np

# Use a single complex dtype for numpy everywhere.
DTYPE = np.complex128

INV_SQRT2 = 1.0 / np.sqrt(2.0)
H = INV_SQRT2 * np.array([[1, 1], [1, -1]], dtype=DTYPE)

# LAMBDA_PI is the base rotation angle realized by the H/T building blocks:
# cos(LAMBDA_PI) = cos^2(pi/8) = (1 + 1/sqrt2)/2. Because LAMBDA_PI / (2 pi) is
# irrational, the multiples {k * LAMBDA_PI mod 2 pi} densely fill [0, 2 pi).
LAMBDA_PI = np.arccos((1.0 + INV_SQRT2) / 2.0)
TWO_PI = 2.0 * np.pi


class Bloch:
    """Axis-angle (Bloch) form of a 2x2 unitary G:

        G = e^{i alpha} (cos(theta/2) I - i sin(theta/2) (n . sigma))

    i.e. a global phase e^{i alpha} times a rotation by angle `theta` about the
    Bloch-sphere axis `n`. Here (n . sigma) = n_x X + n_y Y + n_z Z.
    """

    alpha: float  # global phase
    n: np.ndarray  # unit rotation axis, shape (3,): [n_x, n_y, n_z]
    theta: float  # rotation angle


def to_bloch(g: np.ndarray) -> Bloch:
    """Recover the Bloch form (alpha, n, theta) of a 2x2 unitary `g`."""
    alpha = 0.5 * np.angle(np.linalg.det(g))
    u = np.exp(-1j* alpha)*g
    c = np.real(np.trace(u)/2)
    theta = 2*np.arccos(c)
    X = np.array([[0,1],[1,0]],dtype=DTYPE)
    Y = np.array([[0,-1j],[1j,0]],dtype=DTYPE)
    Z = np.array([[1,0],[0,-1]],dtype=DTYPE)
    s = np.sin(theta/2)
    nx = np.real((-1j/2) * np.trace(X @ u)) / s
    ny = np.real((-1j/2) * np.trace(Y @ u)) /s
    nz = np.real((-1j/2) * np.trace(Z @ u)) /s
    n = np.array([nx, ny, nz])
    n /= np.linalg.norm(n)
    return Bloch(alpha, n, theta)


# n1, n2 are two orthogonal Bloch-sphere axes (n1 . n2 == 0)
# TODO: fill in the two orthogonal rotation axes (each a length-3
# unit vector [x, y, z])
c = np.sqrt(2) / np.tan(np.pi / 8)
n1 = np.array([-c / np.sqrt(2),1, c / np.sqrt(2)], dtype=float)
n1 /= np.linalg.norm(n1)
n2 = np.array([1 / np.sqrt(2),c,-1 / np.sqrt(2)], dtype=float)
n2 /= np.linalg.norm(n2)

# frame derived from the axes (given)
# take the dot product of the Bloch axis with these
# the minus sign arises from the double cover issue
a1 = -n1
a2 = -n2
a3 = np.cross(a1, a2)


def n1n2n1_angles(b: Bloch) -> tuple[float, float, float, float]:
    """Factor the rotation part of a unitary (given as its Bloch form `b`) as
        u = e^{i global_phase} * Rn1(alpha) * Rn2(beta) * Rn1(gamma)

    where Ra(angle) is a rotation by `angle` about axis a, and {a1, a2, a3} is
    the orthonormal frame defined above. Returns (alpha, beta, gamma, global_phase).
    """
    # TODO(student): implement using the steps above.
    v1 = np.dot(b.n,a1)*np.sin(b.theta)
    v2 = np.dot(b.n,a2)*np.sin(b.theta)
    v3 = np.dot(b.n,a3)*np.sin(b.theta)
    beta = np.arctan2(np.sqrt(v2*v2 + v3*v3),np.cos(b.theta))
    sum = np.arctan2(v1,np.cos(b.theta))
    diff = np.arctan2(v3,v2)
    alpha = (sum - diff)/2
    gamma = (sum + diff)/2
    return alpha, beta, gamma, b.alpha


def approx_angle_with_tolerance(angle: float, tolerance: float) -> int:
    """Find an integer multiple k such that
        (k * LAMBDA_PI) mod 2*pi  ~=  angle   (within `tolerance`)
    Since LAMBDA_PI / (2 pi) is irrational, such a k always exists; search
    k = 1, 2, 3, ... and return the first one whose wrapped multiple lands within
    `tolerance` of `angle` (compare both as angles in [0, 2 pi)).

    Hint:
      * wrap an angle into [0, 2 pi)
      * the angular distance between two wrapped angles a, b is
        min(|a - b|, TWO_PI - |a - b|) (so 0.01 and 2*pi - 0.01 count as close).
    """
    # TODO(student): implement using the hint above.
    target = angle % TWO_PI
    k = 1
    while True:
        current = (k * LAMBDA_PI) % TWO_PI
        diff = abs(current - target)
        diff = min(diff, TWO_PI - diff)
        if diff <= tolerance:
            return k
        k += 1



def decompose_2x2(u: np.ndarray, tolerance: float) -> tuple[int, int, int]:
    """Approximate a 2x2 unitary `u` as a product of powers of M1 and M2:

        u  ~=  M1^k * M2^l * M1^m     (up to a global phase)

    where M1 is a rotation about axis a1 and M2 a rotation about axis a2, each by
    the base angle realized by the H/T building blocks. Returns the powers
    (k, l, m).

    Steps (combine the two functions above):

      1. Get the Bloch form of u (to_bloch), then factor its rotation into the
         three frame angles with n1n2n1_angles:
             alpha, beta, gamma, _global_phase = n1n2n1_angles(to_bloch(u))
         alpha and gamma are rotations about a1 (realized by powers of M1);
         beta is a rotation about a2 (realized by powers of M2).

      2. Convert each angle to an integer power with approx_angle_with_tolerance:
             k = approx_angle_with_tolerance(alpha, tolerance)   # power of M1
             l = approx_angle_with_tolerance(beta,  tolerance)   # power of M2
             m = approx_angle_with_tolerance(gamma, tolerance)   # power of M1
         (Mind the relationship between a target rotation angle and the base
         angle each application of M1/M2 adds.)

      3. Return (k, l, m).
    """
    # TODO(student): implement using the steps above.
    b = to_bloch(u)
    alpha, beta, gamma, _ = n1n2n1_angles(b)
    k = approx_angle_with_tolerance(alpha,tolerance)
    l = approx_angle_with_tolerance(beta,tolerance)
    m = approx_angle_with_tolerance(gamma,tolerance)
    return k, l, m

def from_axis_angle(b: Bloch) -> np.ndarray:
    """Build a 2x2 unitary from its Bloch form: a global phase times a rotation
    by angle b.theta about axis b.n (inverse of to_bloch).

        G = e^{i b.alpha} (cos(b.theta/2) I - i sin(b.theta/2) (b.n . sigma))

    where (b.n . sigma) = n_x X + n_y Y + n_z Z. Assumes b.n is a unit vector.
    """
    # TODO: implement using the formula above.
    I = np.eye(2,dtype=DTYPE)
    X = np.array([[0,1],[1,0]],dtype=DTYPE)
    Y = np.array([[0,-1j],[1j,0]],dtype=DTYPE)
    Z = np.array([[1,0],[0,-1]],dtype=DTYPE)
    sigma = b.n[0]*X + b.n[1]*Y + b.n[2]*Z
    c = np.cos(b.theta/2)
    s = np.sin(b.theta/2)
    return np.exp(1j*b.alpha)*(c*I - 1j*s *sigma)

def Rz(theta: float) -> np.ndarray:
    """Rotation about the z axis (no global phase):

    Rz(theta) = diag(e^{-i theta/2}, e^{i theta/2}).
    """
    # TODO: implement (hint: from_axis_angle about axis [0, 0, 1]).
    return from_axis_angle(Bloch(alpha=0.0,n=np.array([0.0, 0.0, 1.0]),theta=theta))

def Ry(theta: float) -> np.ndarray:
    """Rotation about the y axis (no global phase):

    Ry(theta) = [[cos(theta/2), -sin(theta/2)], [sin(theta/2), cos(theta/2)]].
    """
    # TODO: implement (hint: from_axis_angle about axis [0, 1, 0]).
    return from_axis_angle(Bloch(alpha=0.0,n=np.array([0.0,1.0,0.0]),theta=theta))

def euler_angles_zyz(u: np.ndarray) -> tuple[float, float, float, float]:
    """ZYZ Euler decomposition of a 2x2 unitary: angles (alpha, beta, gamma, delta)
    with

        u = e^{i alpha} Rz(beta) Ry(gamma) Rz(delta).

    alpha is the global phase (arg(det u)/2); the rest come from S = e^{-i alpha} u
    in SU(2), where s00 = cos(gamma/2) e^{-i(beta+delta)/2} and
    s10 = sin(gamma/2) e^{i(beta-delta)/2}. When gamma = 0 (s10 = 0), beta/delta are
    split arbitrarily (gimbal lock); the identity still holds.
    """
    # TODO: implement using the relations above.
    alpha = 0.5*np.angle(np.linalg.det(u))
    S = np.exp(-1j*alpha)*u
    gamma = 2 *np.arctan2(abs(S[1,0]), abs(S[0,0]))
    if abs(S[1,0]) <pow(10,-12):
        beta = -2*np.angle(S[0,0])
        delta = 0.0
    else:
        sum = -2*np.angle(S[0,0]) 
        diff = 2*np.angle(S[1,0])     
        beta = (sum+diff) /2
        delta = (sum-diff)/ 2
    return alpha,beta,gamma,delta

def unitary2_sqrt(u: np.ndarray) -> np.ndarray:
    """Principal square root: a 2x2 unitary V with V @ V == u, phase included.
    Take the Bloch form of u and halve both alpha and theta (same axis); squaring
    back doubles them, reproducing u exactly.
    """
    # TODO: implement (hint: to_bloch, halve alpha and theta, from_axis_angle).
    b = to_bloch(u)
    return from_axis_angle(Bloch(alpha=b.alpha/2,n=b.n,theta=b.theta/2) )


# ---------------------------------------------------------------------------
# H/T word machinery for approximating a 2x2 unitary in {H, T} (see cpp/src/HT.h).
#
# M1, M2 are short H/T words that realize rotations by THETA_M = 2*LAMBDA_PI about
# the axes a1, a2. A word is a flat string of 'H'/'T' characters, read left-to-right
# as a matrix product (leftmost char = leftmost/outermost factor).
# ---------------------------------------------------------------------------

# alternating (T-power, H-power, ...) exponents, starting with T
M1_WORD = [7, 1, 1, 1]
M2_WORD = [2, 1, 1, 1, 6, 1, 7, 1, 5, 1, 1, 1, 2, 1, 1, 1, 2, 1, 7, 1, 6]


def expand_word(word: list[int]) -> str:
    """Flatten an alternating (T-power, H-power, ...) exponent list into a literal
    string of 'H'/'T' gates (left-to-right). Even indices are T, odd indices are H.
    """
    # TODO: implement.
    ans = ""
    for i, p in enumerate(word):
        gate = "T" if i%2 == 0 else "H"
        ans+=gate*p
    return ans

# flat H/T strings for the two building-block words (computed once expand_word works)
M1_STR = expand_word(M1_WORD)
M2_STR = expand_word(M2_WORD)


def gates_to_unitary(gates: str) -> np.ndarray:
    """The 2x2 unitary of a flat H/T gate string (left-to-right product)."""
    # TODO: implement (multiply H / T for each char, starting from I).
    T=np.array([[1,0],[0,np.exp(1j*np.pi/4)]],dtype=DTYPE)
    U=np.eye(2,dtype=DTYPE)
    for g in gates:
        if g=="H":
            U=np.dot(U,H)
        elif g=="T":
            U=np.dot(U,T)
    return U


def invert_gates(gates: str) -> str:
    """Inverse of a flat H/T word: reverse the gate order and invert each gate.
    H^-1 = H; the {H, T} basis has no T-dagger, so T^-1 must be spelled as T^7.
    """
    # TODO: implement.
    ans=""
    for g in reversed(gates):
        if g=="H":
            ans+="H"
        elif g=="T":
            ans+="TTTTTTT"
    return ans



def power_gates(base: str, k: int) -> str:
    """The k-th power of a flat H/T word: base repeated k times. Negative k uses the
    inverse word (invert_gates).
    """
    # TODO: implement.
    if k==0:
        return ""
    if k>0:
        return base*k
    return invert_gates(base)*(-k)


def approximate_in_ht(u: np.ndarray, error: float) -> str:
    """Approximate a 2x2 unitary `u` by a flat H/T word (up to global phase) to the
    angular tolerance `error` (smaller -> longer, more accurate).

    Use decompose_2x2 to get the powers (k, l, m) with u ~= M1^k M2^l M1^m, then
    assemble the word:

        power_gates(M1_STR, k) + power_gates(M2_STR, l) + power_gates(M1_STR, m).
    """
    # TODO: implement using decompose_2x2 and power_gates.
    k,l,m=decompose_2x2(u,error)
    return power_gates(M1_STR,k)+power_gates(M2_STR,l)+power_gates(M1_STR,m)
