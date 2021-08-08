import pytest

import numpy as np

import BondGraphTools as bgt
from BondGraphTools.reaction_builder import Reaction_Network


def four_site_phos_model():
    rn = Reaction_Network(name="4-site sequential phosphorylation")
    rn.add_reaction("E + S0 + ATP = ES0", name="e1a")
    rn.add_reaction("ES0 = E + S1 + ADP", name="e1b")
    rn.add_reaction("F + S1 = FS1", name="f1a")
    rn.add_reaction("FS1 = F + S0 + Pi", name="f1b")
    rn.add_reaction("E + S1 + ATP = ES1", name="e2a")
    rn.add_reaction("ES1 = E + S2 + ADP", name="e2b")
    rn.add_reaction("F + S2 = FS2", name="f2a")
    rn.add_reaction("FS2 = F + S1 + Pi", name="f2b")
    rn.add_reaction("E + S2 + ATP = ES2", name="e3a")
    rn.add_reaction("ES2 = E + S3 + ADP", name="e3b")
    rn.add_reaction("F + S3 = FS3", name="f3a")
    rn.add_reaction("FS3 = F + S2 + Pi", name="f3b")
    rn.add_reaction("E + S3 + ATP = ES3", name="e4a")
    rn.add_reaction("ES3 = E + S4 + ADP", name="e4b")
    rn.add_reaction("F + S4 = FS4", name="f4a")
    rn.add_reaction("FS4 = F + S3 + Pi", name="f4b")

    rn.add_chemostat("ATP")
    rn.add_chemostat("ADP")
    rn.add_chemostat("Pi")

    model = rn.as_network_model(normalised=True)
    return model


def set_4site_model_parameters(system):
    # Parameterise chemical potentials of ATP, ADP and Pi to be approximately
    # consistent with physiological free energy of around -50 to -60 kJ/mol
    ATP_potential = 50000 / 8.314 / 310  # ΔG_ATP/(RT) (dimensionless)
    ADP_potential = 0.0
    Pi_potential = 0.0

    # Assume that the affinities of K_E and K_F are 1
    K_E = 1.0
    K_F = 1.0

    args = (ATP_potential, ADP_potential, Pi_potential, K_E, K_F)

    gamma = 34700.0

    K_Sn, r_Ea, r_Eb, K_ESn, r_Fa, r_Fb, K_FSn = \
        foursite_energetic_params(gamma, *args)

    # Set parameters for bond graph model
    (system / "SS:ATP").set_param('e', ATP_potential)
    (system / "SS:ADP").set_param('e', ADP_potential)
    (system / "SS:Pi").set_param('e', Pi_potential)
    (system / "C:E").set_param('k', K_E)
    (system / "C:F").set_param('k', K_F)

    for i in range(5):
        (system / f"C:S{i}").set_param('k', K_Sn[i])

    for i in range(4):
        (system / f"C:ES{i}").set_param('k', K_ESn[i])
        (system / f"R:e{i+1}a").set_param('r', r_Ea[i])
        (system / f"R:e{i+1}b").set_param('r', r_Eb[i])
        (system / f"C:FS{i+1}").set_param('k', K_FSn[i])
        (system / f"R:f{i+1}a").set_param('r', r_Fa[i])
        (system / f"R:f{i+1}b").set_param('r', r_Fb[i])


def foursite_energetic_params(
        gamma, ATP_potential, ADP_potential, Pi_potential, K_E, K_F):

    # Kinetic constants
    aE = [0.00812, 0.102, 0.00812, 0.102]  # unit nM/sec
    bE = [0.016, 0.204, 0.016, 0.204]  # unit 1/sec
    cE = [0.100, 10.00, 0.100, 10.000]  # unit 1/sec
    aF = [0.112, 0.00264, 0.651, 0.00285]  # unit nM/sec
    bF = [0.224, 0.005, 1.30, 0.006]  # unit 1/sec
    cF = [11.0, 0.017, 63.9, 0.136]  # unit 1/sec

    # Calculate affinities of protein substrates.
    # Assume that they follow the equation
    # K_Sn = γ^n
    K_Sn = [gamma**i for i in range(5)]

    # Calculate affinities of the rest of the energetic parameters
    r_Ea = [aE[i] / (K_Sn[i] * K_E * np.exp(ATP_potential)) for i in range(4)]
    K_ESn = [b / ra for b, ra in zip(bE, r_Ea)]
    r_Eb = [c / K_ES for c, K_ES in zip(cE, K_ESn)]
    r_Fa = [aF[i] / (K_Sn[i + 1] * K_F) for i in range(4)]
    K_FSn = [b / ra for b, ra in zip(bF, r_Fa)]
    r_Fb = [c / K_FS for c, K_FS in zip(cF, K_FSn)]

    return K_Sn, r_Ea, r_Eb, K_ESn, r_Fa, r_Fb, K_FSn


def test_model_vars():
    model = four_site_phos_model()
    assert len(model.state_vars) == 15
    assert len(model.control_vars) == 34


@ pytest.mark.slow
def test_simulation():
    model = four_site_phos_model()
    set_4site_model_parameters(model)

    # Initial concentrations (unit nM)
    S_tot = 1e4
    E_tot = 2.8e3
    F_tot = 2.8e3

    species = [comp.name for comp, q in model.state_vars.values()]
    idx_E = species.index("E")
    idx_F = species.index("F")
    idx_S0 = species.index("S0")
    idx_S4 = species.index("S4")

    alpha_vals = [0.0, 0.5, 0.9]
    expected_S4 = [6077.03288712802, 2590.0823948000684, 6.480921602618817e-07]

    for alpha, S4_end in zip(alpha_vals, expected_S4):
        x0 = np.zeros(len(model.state_vars))
        x0[idx_E] = E_tot
        x0[idx_F] = F_tot
        x0[idx_S0] = alpha * S_tot
        x0[idx_S4] = (1 - alpha) * S_tot

        # Use different time steps for small and large time
        tspan = (0.0, 1000.0)
        t, x = bgt.simulate(model, tspan, x0, dt=1)
        assert x[-1][idx_S4] == pytest.approx(S4_end, abs=1e-6)
