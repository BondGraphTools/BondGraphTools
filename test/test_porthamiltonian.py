import BondGraphTools as bgt
from BondGraphTools.port_hamiltonian import PortHamiltonian
import sympy as sp


def test_hamiltonian():
    hamiltonian = "x^2/(2*C) + x^4/(4*D)"

    relations, state_vars, params, ports = PortHamiltonian._generate_relations(
        hamiltonian
    )

    assert state_vars == {"q_0": "x"}
    assert params == {"C": None, "D":None}
    assert ports == {0:None}
    assert relations == ["-e_0 + q_0**3/D + q_0/C", "dq_0 - f_0"]


def test_create_PH():
    hamiltonian = "x^2/(2*C) + x^4/(4*D)"
    ph = bgt.new(component="PH", value=hamiltonian)
    assert ph.params == {"C": None, "D":None}

    assert ph.state_vars == {"q_0": "x"}

    assert sp.sympify("-e_0 + q_0**3/D + q_0/C") in ph.constitutive_relations
    assert sp.sympify("dq_0 - f_0") in ph.constitutive_relations

    assert ph.hamiltonian == hamiltonian
    port_list = list(ph.ports)
    assert len(port_list) == 1
    assert port_list[0] == (ph, 0)


def test_create_PH_2():
    hamiltonian = "x^2/(2*C + 2*y)"
    build_args = {
        "hamiltonian": hamiltonian,
    }
    ph = bgt.new(component="PH", value=build_args)

    assert ph.params == {"C": None}

    assert set(ph.state_vars.values()) == {"x", "y"}

    for rel in ("-e_0 + q_0/(C + q_1)",
                "dq_0 - f_0",
                "-e_1 - q_0**2/(2*(C + q_1)**2)",
                "dq_1 - f_1"):
        assert sp.sympify(rel) in ph.constitutive_relations


    assert ph.hamiltonian == hamiltonian

    port_list = list(ph.ports)
    assert len(port_list) == 2
    assert (ph, 0) in port_list
    assert (ph, 1) in port_list


def test_create_PH_parameters():
    hamiltonian = "x^2/(2*C) + x^4/(4*D)"
    p_1 = 1
    p_2 = sp.S("k")
    build_args = {
        "hamiltonian": hamiltonian,
        "params":{"C":p_1, "D":p_2}
    }

    ph = bgt.new(component="PH", value=build_args)
    assert ph.params == {"C": p_1, "D": p_2}

    assert ph.state_vars == {"q_0": "x"}

    assert sp.sympify("-e_0 + q_0**3/k + q_0") in ph.constitutive_relations
    assert sp.sympify("dq_0 - f_0") in ph.constitutive_relations
    assert ph.hamiltonian == hamiltonian


def test_create_duffing_eqn():

    model = bgt.new(name="Duffing")
    junct = bgt.new("0")
    nlin_cap = bgt.new("PH", value="x^2/2 + x^4/4")
    inductor = bgt.new("I", value=1)
    source = bgt.new("Sf")
    bgt.add(model, nlin_cap, inductor, source, junct)
    bgt.connect(junct, nlin_cap)
    bgt.connect(junct, inductor)
    bgt.connect(source, junct)

    rel = model.constitutive_relations

    assert sp.sympify("dx_0 + x_1 - u_0") in rel

    assert sp.sympify("dx_1 - x_0**3 - x_0") in rel
