from sympy import SparseMatrix, symbols, Pow, Matrix
from collections import defaultdict
import logging
from math import log
from .actions import connect, disconnect, new
logger = logging.getLogger(__name__)

# Gas constant, J/Mo/K
R = 8.3144598

# Component Library
LIBRARY = "BioChem"


class Reaction_Network(object):
    def __init__(self, reactions=None, name=None, temperature=300, volume=1):
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
        self._T = temperature
        self._V = volume
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
    def species(self):
        return list(self._species.keys())

    @property
    def stoichiometry(self):

        return self.reverse_stoichiometry - self.forward_stoichiometry

    @property
    def reverse_stoichiometry(self):
        matrix = SparseMatrix(len(self._species), len(self._reactions),{})
        for col, (_, forward_species, _, _) in enumerate(
                self._reactions.values()):
            for species, qty in forward_species.items():
                matrix[(self.species.index(species), col)] = qty

        return matrix

    @property
    def forward_stoichiometry(self):
        matrix = SparseMatrix(len(self._species), len(self._reactions),{})
        for col, (back_species,_, _, _) in enumerate(
                self._reactions.values()):
            for species, qty in back_species.items():
                matrix[(self.species.index(species), col)] = qty

        return matrix

    def fluxes(self):
        coords = symbols(",".join([f"x_{{{s}}}" for s in self.species]))
        fluxes = []
        for col, (back_species, forward_species, _, _) in enumerate(
                self._reactions.values()):
            prod = 1
            subst = 1
            for s, qty in forward_species.items():
                subst *= Pow(symbols(f"k_{{{s}}}") *
                             symbols(f"x_{{{s}}}"), qty)

            for s, qty in back_species.items():
                prod *= Pow(symbols(f"k_{{{s}}}") *
                            symbols(f"x_{{{s}}}"), qty)
            eqn = prod - subst
            fluxes.append(eqn)
        matrix = Matrix(len(fluxes),1, fluxes )
        return matrix, coords

    def as_network_model(self, normalised=False):

        system = new(name=self.name)
        species_anchor = self._build_species(system, normalised)
        self._build_reactions(system, species_anchor, normalised)

        return system

    def _build_reactions(self, system, species_anchors, normalised):

        if normalised:
            param_dict = {"R": 1, "T": 1}
        else:
            param_dict = {"R": R, "T": self._T}

        for reaction_name, (bck_sto, fwd_sto, _, _) in self._reactions.items():
            reaction = new("Re", library=LIBRARY,
                           name=reaction_name, value=param_dict)
            system.add(reaction)
            fwd_name = "".join(
                list(fwd_sto.keys())
            )
            bck_name = "".join(
                list(bck_sto.keys())
            )
            if len(bck_sto) == 1 and list(bck_sto.values())[0] == 1:
                species = list(bck_sto.keys()).pop()
                connect(reaction, species_anchors[species])
            else:
                reverse_complex = new("1", name=bck_name)
                system.add(reverse_complex)
                connect(reaction, reverse_complex.inverting)
                self._connect_complex(
                    system, species_anchors, reverse_complex, bck_sto
                )
            if len(fwd_sto) == 1 and list(fwd_sto.values())[0] == 1:
                species = list(fwd_sto.keys()).pop()
                connect(reaction, species_anchors[species])
            else:
                forward_complex = new("1", name=fwd_name)
                system.add(forward_complex)

                connect(reaction, forward_complex.inverting)

                self._connect_complex(
                    system, species_anchors, forward_complex, fwd_sto
                )

    def _connect_complex(self, system, species_anchors, junct, stoichiometry):
        for i, (species, qty) in enumerate(stoichiometry.items()):

            if qty == 1:
                connect(junct.non_inverting, species_anchors[species])
            else:
                tf = new("TR", value=qty)
                system.add(tf)
                connect((tf, 1), junct.non_inverting)
                connect(species_anchors, (tf,0))

    def _build_species(self, system, normalised):
        if normalised:
            param_dict = {"R": 1, "T": 1}
        else:
            param_dict = {"R": R, "T": self._T}

        species_anchors = {}
        for species, n_reactions in self._species.items():
            # edit: create new component for chemostat
            if species in self._chemostats:
                this_species = new(
                        "Se", name=species, value=self._chemostats[species]
                )
                n_reactions = n_reactions-1
            else:
                this_species = new(
                        "Ce", library=LIBRARY, name=species, value=param_dict
                )
            system.add(this_species)

            if n_reactions == 1:
                species_anchors[species] = this_species
            else:
                anchor = new("0", name=species)
                system.add(anchor)
                connect(this_species, anchor)
                species_anchors[species] = anchor

            if species in self._flowstats:
                flowstat = new(
                    "Sf", value=self._flowstats[species], name=species
                )
                system.add(flowstat)
                connect(flowstat, species_anchors[species])


#            if species in self._chemostats:
#                chemostat = new("Se", value=0, name=species)
#                connect(chemostat, species_anchors[species])
#                this_species.initial_values['p_0'] = self._chemostats[species]

        return species_anchors

    def add_reaction(self, reaction,
                     forward_rates=None, reverse_rates=None, name=""):

        reaction_step = []
        remaining_reactions = reaction

        if not name or name in self._reactions:
            n = 1
            while "r_{{{base}}}".format(base=str(n)+';0') in self._reactions:
                n += 1
            base = str(n) + ';'
            idx = "r_{{{base}}}"
        else:
            idx = "{name}".format(name=name)
            base = ""

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

            self._reactions[idx.format(base=base+str(i))] = (reaction_step[i],
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
        if species not in self._flowstats:
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

