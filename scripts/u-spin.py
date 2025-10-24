"""Calculate the U-spin limit for ACP(KK) and ACP(pipi), taking into account the fact that their branching fractions
differ.

    ACP(KK) / ACP(pipi) ~ - A(pipi) / A(KK)

where A ~ sqrt(BF / PHSP)
"""


def kallen(a, b, c):
    return a**2 + b**2 + c**2 - 2 * a * b - 2 * a * c - 2 * b * c


def momentum_com_2body(M, m1, m2):
    """Return the momentum in the CoM frame for a 2-body decay of a particle of mass M."""
    return (kallen(M**2, m1**2, m2**2)) ** 0.5 / (2 * M)


if __name__ == "__main__":
    bf_pp = 1.453e-3
    bf_kk = 4.08e-3

    m_d0 = 1864.84
    m_pi = 139.57039
    m_k = 493.677

    p_pp = momentum_com_2body(m_d0, m_pi, m_pi)
    p_kk = momentum_com_2body(m_d0, m_k, m_k)

    # This is the ratio of "tree" (read CKM dominant) amplitudes: amplitude_ratio = tree_kk / tree_pipi
    # Then aCP(kk) / aCP(pipi) * amplitude_ratio = -1 in the U-spin limit
    amplitude_ratio = ((bf_kk / p_kk) / (bf_pp / p_pp)) ** 0.5

    print(f"amplitude_ratio A_KK / A_pipi = {amplitude_ratio:.4f}")
