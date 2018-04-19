import BondGraphTools as bgt
from BondGraphTools.sim_tools import simulate

def test_c_sim():

    c = bgt.new("C", value=0.001)
    r = bgt.new("R", value=1000)
    Se = bgt.new("Se")
    kcl = bgt.new("0")

    bg = c + r + Se + kcl

    bg.connect([(c, kcl),
                (r, kcl),
                (Se, kcl)])

    [t, x] = simulate(
        bg,
        intial=[0],
        params=["exp(-t)"],
        observable=[Se.ports, c.state_vars]
    )
