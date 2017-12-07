from .reactionparser import create_reaction, _split_reactants, \
    join_bond_graphs, reaction_to_bond_graph



def test_parser():
    process = "A + 2*B = C + 3*D"
    nodes, edges = create_reaction(process).pop()

def test_join():
    process_1 = "A + 2*B = C + 3*D"
    process_2 = "A = C + E"

    nodes_1, edges_1 = create_reaction(process_1).pop()
    nodes_2, edges_2 = create_reaction(process_2).pop()

    nodes_3, edges_3 = join_bond_graphs(
        (nodes_1, edges_1), (nodes_2, edges_2)
    )

    for k, v in nodes_3.items():
        if k in nodes_1:
            assert v == nodes_1[k]

        if k in nodes_2:
            assert v == nodes_2[k]

        if k not in nodes_1 and k not in nodes_2:
            assert False, "Key {k} not found".format(k=k)


def test_split_reactants():

    reactants = "A + 2*B"
    res = _split_reactants(reactants)

    assert res == {"A": 1, "B": 2}