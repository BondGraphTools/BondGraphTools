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
        cat_reaction = bgt.new("Re", name="Re", library="BioChem", value={'r': None, 'R': 1, 'T': 1})

        # We choose 'k' to be 1 for demonstration.
        enzyme = bgt.new("Ce", name="E", library="BioChem", value={'k': 1, 'R': 1, 'T': 1})

        # Substrate + Enzyme flux conservation law
        SE = bgt.new('1')
        # Product + Enzyme flux conservation law
        PE = bgt.new('1')

        # Conservation of enzyme law.
        law_E = bgt.new("0")

        bgt.add(cat_model, substrate, product, enzyme, SE, PE, law_E, cat_reaction)

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