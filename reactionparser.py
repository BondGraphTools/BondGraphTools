

def reaction_to_bond_graph(forward_affinities, reverse_affinities, reaction_id):

    tf_count = 0
    reaction = "r{index}".format(index=reaction_id)
    reaction_in = "1: Re{index}-In".format(index=reaction_id)
    reaction_out = "1: Re{index}-Out".format(index=reaction_id)
    nodes = {
        reaction: "Re: " + reaction,
        reaction_in: "1",
        reaction_out: "1"
    }

    edges = [(reaction_in, reaction), (reaction, reaction_out)]

    for target, affinities in [(reaction_in, forward_affinities),
                               (reaction_out, reverse_affinities)]:
        for species, stoich_coeff in affinities.items():

            node = "0: {species}".format(species=species)

            if node not in nodes:
                source = "C: {species}".format(species=species)
                nodes[source] = source
                nodes[node] = "0"
                edges.append(
                    (node, source)
                )

            if stoich_coeff != 1:
                tf_count += 1
                transformer = "TF: {reaction_id},{tf_count}".format(
                    tf_count=tf_count,
                    reaction_id=reaction_id
                )
                nodes[transformer] = "TF: {c}".format(c=stoich_coeff)
                edges.append(
                    (node, transformer)
                )
                edges.append(
                    (transformer, target)
                )
            else:
                edges.append(
                    (node, target)
                )
    return nodes, edges


def join_bond_graphs(*args):
    """
    Merges parallel bond graphs together.

    Args:
        args (iterable): An iterable of (nodes, edges)

    Nodes and edges are of the form ::
        nodes = {"node_id": "display string"}
        edges =[("base_node_id", "target_node_id")]

    """

    nodes = {}
    edges = []

    for (arg_nodes, arg_edges) in args:
        nodes.update(arg_nodes)
        edges += arg_edges

    return nodes, edges


def create_reaction(reaction, forward_rates=None, reverse_rates=None):
    """
    Creates an instance of :obj:`Reaction` corresponding to the argument
    `reaction`.
    The reaction string is expected to be in the form
        "A + B = C + D"

    Further, it is assumed that each species is represented by a sequence of
    alphanumerics, with the first letter being uppercase.

    If forwards and backwards rates are specified; it is expected that the
     number of rates are the same as the number of reactions.

    Args:
        reaction (str): The chemical formula for the reaction.
        forward_rates (float or list(float)): The rate(s) at which
         the reaction(s) proceed(s).
        reverse_rates (float or list(float)): The rate(s) at which the
         reaction(s) reverse(s).

    Returns:
        list(:obj:`Reaction`)

    Raises:

    """

    reaction_step = []
    remaining_reactions = reaction
    reaction_objs =[]
    while remaining_reactions:
        reactants, _, remaining_reactions = remaining_reactions.partition("=")
        stoichiometrics = _split_reactants(reactants)
        reaction_step.append(stoichiometrics)

    for i in range(len(reaction_step) - 1):
        try:
            f_rate = forward_rates[i]
        except TypeError:
            f_rate = forward_rates
        try:
            r_rate = reverse_rates[i]
        except TypeError:
            r_rate = reverse_rates

        reaction_objs.append(
            reaction_to_bond_graph(forward_affinities=reaction_step[i],
                                   reverse_affinities=reaction_step[i+1],
                                   reaction_id=i)
        )
    return reaction_objs


def _split_reactants(reactants):
    reactants = reactants.replace(" ", "").split("+")

    stoiciometrics = dict()

    for reactant in reactants:
        try:
            coeff, prod = reactant.split("*")
            coeff = int(coeff)
        except ValueError:
            prod = reactant
            coeff = 1

        stoiciometrics[prod] = coeff

    return stoiciometrics

