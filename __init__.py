

def rlc_example():
    import model as mo


    bg = mo.BondGraph()

    c = bg.add_component(
        component="C",
        name="C1",
        pos=(0, 1)
    )

    r = bg.add_component(
        component="R",
        name="R1",
        pos=(1, 1)
    )

    i = bg.add_component(
        component="I",
        name="I1",
        pos=(-1, 1)

    )

    kvl = bg.add_component(component="0", name="KVL1", pos=(0, 0))

    bg.add_bond("KVL1", "C1")
    bg.add_bond("KVL1", "R1")
    bg.add_bond("KVL1", "I1")



