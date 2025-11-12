from helpers.ndqd_amplitudes import rMM, tMM, mMM, aMM, rOrth
import qutip as qt


def get_reflection_amplitudes(
    sh: int,
    N: int,
    g: float,
    fP: float,
    fA: float,
    fC: float,
    mm_fc: complex,
    mm_fr: complex,
    k: float,
    kr: float,
    kt: float,
    km: float,
    gamma: float,
    alpha: float,
):
    r_pi = qt.coherent(
        sh,
        rMM(
            N=N,
            g=g,
            fP=fP,
            fA=fA,
            fC=fC,
            mmFC=mm_fc,
            mmFR=mm_fr,
            k=k,
            kr=kr,
            gamma=gamma,
        )
        * alpha,
    )
    print(r_pi)

    t_pi = qt.coherent(
        sh,
        tMM(
            N=N,
            g=g,
            fP=0,
            fA=fA,
            fC=fC,
            mmFC=mm_fc,
            k=k,
            kr=kr,
            kt=kt,
            gamma=gamma,
        )
        * alpha,
    )

    m_pi = qt.coherent(
        sh,
        mMM(
            N=0,
            g=g,
            fP=0,
            fA=fA,
            fC=fC,
            mmFC=mm_fc,
            k=k,
            kr=kr,
            km=km,
            gamma=gamma,
        )
        * alpha,
    )

    a_pi = qt.coherent(
        sh,
        aMM(
            N=0,
            g=g,
            fP=0,
            fA=fA,
            fC=fC,
            mmFC=mm_fc,
            k=k,
            kr=kr,
            gamma=gamma,
        )
        * alpha,
    )

    rO_pi = qt.coherent(
        sh,
        rOrth(
            N=0,
            g=g,
            fP=0,
            fA=fA,
            fC=fC,
            mmFC=mm_fc,
            mmFR=mm_fr,
            k=k,
            kr=kr,
            gamma=gamma,
        )
        * alpha,
    )

    r_v = qt.coherent(
        sh,
        rMM(
            N=N,
            g=g,
            fP=fP + 500,
            fA=fA,
            fC=fC,
            mmFC=mm_fc,
            mmFR=mm_fr,
            k=k,
            kr=kr,
            gamma=gamma,
        )
        * alpha,
    )

    t_v = qt.coherent(
        sh,
        tMM(
            N=N,
            g=g,
            fP=fP + 500,
            fA=fA,
            fC=fC,
            mmFC=mm_fc,
            k=k,
            kr=kr,
            kt=kt,
            gamma=gamma,
        )
        * alpha,
    )

    m_v = qt.coherent(
        sh,
        mMM(
            N=N,
            g=g,
            fP=fP + 500,
            fA=fA,
            fC=fC,
            mmFC=mm_fc,
            k=k,
            kr=kr,
            km=km,
            gamma=gamma,
        )
        * alpha,
    )

    a_v = qt.coherent(
        sh,
        aMM(
            N=N,
            g=g,
            fP=fP + 500,
            fA=fA,
            fC=fC,
            mmFC=mm_fc,
            k=k,
            kr=kr,
            gamma=gamma,
        )
        * alpha,
    )

    rO_v = qt.coherent(
        sh,
        rOrth(
            N=N,
            g=g,
            fP=fP + 500,
            fA=fA,
            fC=fC,
            mmFC=mm_fc,
            mmFR=mm_fr,
            k=k,
            kr=kr,
            gamma=gamma,
        )
        * alpha,
    )

    pi = [r_pi, t_pi, m_pi, a_pi, rO_pi]
    v = [r_v, t_v, m_v, a_v, rO_v]

    return pi, v
