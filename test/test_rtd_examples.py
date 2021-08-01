def test_example_2():
    import BondGraphTools as bgt

    def enzyme_catalysed_reaction(name):
        """
        This function produces a bond graph model of an basic enzyme catalysed
        reaction of the from `S + E  = E + P` where the substrate and product
        are exposed as external ports.

        Args:
            name (str): The name of the enzyme

        Returns:
            `BondGraph`: The resulting model
        """

        cat_model = bgt.new(name=name)

        # Construct the external ports.
        substrate = bgt.new("SS", name="S")
        product = bgt.new("SS", name="P")

        # Here we build the reaction, again with the rate as a control variable.
        # Again, we assume parameterised have be normalised with respect to
        # pressure and temperature.
        cat_reaction = bgt.new(
            "Re", name="Re", library="BioChem", value={
                'r': None, 'R': 1, 'T': 1})

        # We choose 'k' to be 1 for demonstration.
        enzyme = bgt.new(
            "Ce",
            name="E",
            library="BioChem",
            value={
                'k': 1,
                'R': 1,
                'T': 1})

        # Substrate + Enzyme flux conservation law
        SE = bgt.new('1')
        # Product + Enzyme flux conservation law
        PE = bgt.new('1')

        # Conservation of enzyme law.
        law_E = bgt.new("0")

        bgt.add(
            cat_model,
            substrate,
            product,
            enzyme,
            SE,
            PE,
            law_E,
            cat_reaction)

        connections = [
            (substrate, SE),
            (law_E, SE),
            (law_E, enzyme),
            (SE, cat_reaction),
            (cat_reaction, PE),
            (PE, law_E),
            (PE, product)
        ]
        for tail, head in connections:
            bgt.connect(tail, head)

        bgt.expose(substrate, 'S')
        bgt.expose(product, 'P')

        return cat_model

    E1 = enzyme_catalysed_reaction("E1")


def test_init_example():

    import BondGraphTools as bgt

    # Create a new model
    model = bgt.new(name="RLC")

    # Create components
    # 1 Ohm Resistor
    resistor = bgt.new("R", name="R1", value=1.0)

    # 1 Henry Inductor
    inductor = bgt.new("I", name="L1", value=1.0)
    # 1 Farad Capacitor
    capacitor = bgt.new("C", name="C1", value=1.0)
    # Conservation Law
    law = bgt.new("0")  # Common voltage conservation law
    bgt.add(model, resistor, inductor, capacitor, law)
    # Connect the components
    bgt.connect(law, resistor)
    bgt.connect(law, capacitor)
    bgt.connect(law, inductor)

    # produce timeseries data
    t, x = bgt.simulate(model, x0=[1, 1], timespan=[0, 10])


def test_rlc_example():

    import BondGraphTools as bgt
    model = bgt.new(name="RC")

    C = bgt.new("C", value=1)
    R = bgt.new("R", value=1)
    zero_law = bgt.new("0")

    bgt.add(model, R, C, zero_law)

    bgt.connect(R, zero_law)
    bgt.connect(zero_law, C)

    bgt.draw(model)

    timespan = [0, 5]
    x0 = [1]
    t, x = bgt.simulate(model, timespan=timespan, x0=x0)
    import matplotlib.pyplot as plt

    fig = plt.plot(t, x)

    Sf = bgt.new('Sf')
    bgt.add(model, Sf)
    bgt.connect(Sf, zero_law)

    bgt.draw(model)

    timespan = [0, 5]
    x0 = [1]
    t, x = bgt.simulate(model, timespan=timespan,
                        x0=x0, control_vars={'u_0': 2})
    plt.plot(t, x)

    t, x = bgt.simulate(model, timespan=timespan, x0=x0,
                        control_vars={'u_0': 'sin(2*t)'})
    plt.plot(t, x)

    def step_fn(t, x, dx):
        return 1 if t < 1 else 0

    t, x = bgt.simulate(model, timespan=timespan, x0=x0,
                        control_vars={'u_0': step_fn})
    plt.plot(t, x)

    fig = plt.figure()
    for i in range(4):
        func_text = "cos({i}*t)".format(i=i)
        t_i, x_i = bgt.simulate(model, timespan=timespan,
                                x0=x0, control_vars={'u_0': func_text})
        plt.plot(t_i, x_i)
