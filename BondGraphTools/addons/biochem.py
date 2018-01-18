from collections import defaultdict
import logging

from ..core.bondgraph import BondGraph


logger = logging.getLogger(__name__)

## Gas constant, J/Mo/K
R = 8.3144598
LIBRARY = "BioChem"

class Reaction_Network(object):

    def __init__(self, reactions=None, name=None, T=300):
        """
        Args:
            reactions:
            name:
            T: Temperature in Kelvin (Default 300K, or approx 27c)
        """
        self._reactions = {}
        self._flowstats = {}
        self._chemostats = {}
        self._species = defaultdict(int)
        self._bond_graph = None
        if name:
            self.name = name
        else:
            self.name = "New Reaction"
        if isinstance(reactions, str):
            self.add_reaction(reactions)
        elif isinstance(reactions, list):
            for reaction in reactions:
                self.add_reaction(reaction)

    @property
    def bond_graph(self):
        """
        The bond graph representation of this particular reaction network.
        """
        if not self._bond_graph:
            self._build_bond_graph()

        return self._bond_graph

    def _build_bond_graph(self):
        graph = BondGraph(self.name)
        species_base = {}
        for species, ref in self._species.items():
            c = graph.add_component(component="Ce",
                                    library=LIBRARY,
                                    name=species)
            if ref == 1:
                species_base[species] = c
            else:
                z = graph.add_component("0")
                graph.add_bond(z, c)
                species_base[species] = z

        for species, concentration in self._chemostats.items():
            cs = graph.add_component(
                'Se', name=species,
                local_params={"e": concentration}
            )
            graph.add_bond(cs, species_base[species])

        for species, flux in self._flowstats.items():
            fs = graph.add_component(
                "Sf", name=species, local_params={"f": flux}
             )
            graph.add_bond(fs, species_base[species])

        for idx in self._reactions:
            re = graph.add_component("Re", name=idx, library=LIBRARY)
            reactants, products, _, _ = self._reactions[idx]

            if len(reactants) == 1:
                in_pipe = {"to_component":re, "to_port": 0}
            else:
                cplex = graph.add_component("1")
                graph.add_bond(cplex, re, None, 0)
                in_pipe = {"to_component":cplex, "to_port": None}

            if len(products) == 1:
                out_pipe = {"from_component":re, "from_port": 1}
            else:
                cplex = graph.add_component("1")
                graph.add_bond(cplex, re, None, 1)
                out_pipe = {"from_component":cplex, "from_port": None}

            for sp, factor in reactants.items():
                if factor == 1:
                    graph.add_bond(species_base[sp], **in_pipe)
                else:
                    tf = graph.add_component("TF", a=factor)
                    graph.add_bond(species_base[sp], tf, to_port=0)
                    graph.add_bond(tf, from_port=1, **in_pipe)

            for sp, factor in products.items():
                if factor == 1:
                    graph.add_bond(**out_pipe, to_component=species_base[sp])
                else:
                    tf = graph.add_component("TF", a=factor)
                    graph.add_bond(**out_pipe, to_component=tf, to_port=0)
                    graph.add_bond(tf, from_port=1,
                                   to_component=species_base[sp])

        self._bond_graph = graph

    def add_reaction(self, reaction,
                     forward_rates=None, reverse_rates=None, name=""):
        reaction_step = []
        remaining_reactions = reaction

        self._bond_graph = None

        if not name or name in self._reactions:
            n = 1

            while "r{}_1".format(n) in self._reactions:
                n += 1

            idx = "r{}_{{}}".format((n))
        else:
            idx = "{name}_{{}}".format(name=name)

        while remaining_reactions:
            in_react, _, remaining_reactions = remaining_reactions.partition(
                "=")
            reactants = _split_reactants(in_react)

            reaction_step.append(reactants)

        for i in range(len(reaction_step) - 1):
            try:
                f_rate = forward_rates[i]
            except TypeError:
                f_rate = forward_rates
            try:
                r_rate = reverse_rates[i]
            except TypeError:
                r_rate = reverse_rates

            for sp in reaction_step[i]:
                self._species[sp] += 1

            for sp in reaction_step[i+1]:
                self._species[sp] += 1

            self._reactions[idx.format(i)] = (reaction_step[i],
                                              reaction_step[i+1],
                                              f_rate, r_rate)

    def add_chemostat(self, species, concentration=None):
        """
        Add a (generalised) Chemostat to the reaction network.
        This provides a variable source of the particular species so as to
        maintain a particular concentration. This can also act as a I/O source.

        Notes:
             Only one chemostat is available per species; so adding a duplicate
             will overwrite the previous chemostatic concentration (if defined)

        Args:
            species: The name/identifier of the particular chemical species
            concentration: (default None) The fixed concentration. Left as a
             free parameter if None
        """
        if species not in self._chemostats:
            self._bond_graph = None
            self._species[species] += 1

        self._chemostats[species] = concentration

    def add_flowstat(self, species, flux=None):
        """
        Adds a (generalised) Flowstat, which provides a flux of the particular
        species in or out of the reaction network.

        Notes:
            Only one flowstat per species, so adding duplicate flowstats will
            overwrite the previous flowstatic fluxes (if defined)

        See Also:
            :func:`add_chemostat`

        Args:
            species: The name/identifier of the chemical species
            flux: (default None) The rate at which this species is
             added/removed. Left as a free paramter if None

        """
        self._bond_graph = None
        if species not in self._flowstats:
            self._bond_graph = None
            self._species[species] += 1

        self._flowstats[species] = flux

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

