from sympy import SparseMatrix
from collections import defaultdict
import logging

from .base import new

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

        return self.forward_stoichiometry - self.reverse_stoichiometry

    @property
    def forward_stoichiometry(self):
        matrix = SparseMatrix(len(self._species), len(self._reactions),{})
        for col, (_, forward_species, _, _) in enumerate(
                self._reactions.values()):
            for species, qty in forward_species.items():
                matrix[(self.species.index(species), col)] = qty

        return matrix

    @property
    def reverse_stoichiometry(self):
        matrix = SparseMatrix(len(self._species), len(self._reactions),{})
        for col, (back_species,_, _, _) in enumerate(
                self._reactions.values()):
            for species, qty in back_species.items():
                matrix[(self.species.index(species), col)] = qty

        return matrix

    def as_network_model(self, normalised=False):

        system = new(name=self.name)

        if normalised:
            param_dict = {"R": 1, "T": 1}
        else:
            param_dict = {"R": R, "T": self._T}

        species_anchor = self._build_species(system, param_dict)
        self._build_reactions(system, species_anchor, param_dict)

        return system

    def _build_reactions(self, system, species_anchors, param_dict):

        for reaction_name, (bck_sto, fwd_sto, _, _) in self._reactions.items():
            reaction = new("Re", library=LIBRARY,
                           name=reaction_name, value=param_dict)
            fwd_name = "".join(
                list(fwd_sto.keys())
            )
            bck_name = "".join(
                list(bck_sto.keys())
            )
            forward_complex = new("Y", library=LIBRARY, name=fwd_name)
            reverse_complex = new("Y", library=LIBRARY, name=bck_name)
            system.add(reaction, forward_complex, reverse_complex)
            system.connect((reaction, 0), (reverse_complex, 0))
            system.connect((reaction, 1), (forward_complex, 0))

            self._connect_complex(
                system, species_anchors, forward_complex, fwd_sto
            )
            self._connect_complex(
                system, species_anchors, reverse_complex, bck_sto
            )

    @staticmethod
    def _connect_complex(system, species_anchors, complex, stoichiometry):
        for i, (species, qty) in enumerate(stoichiometry.items()):
            port = i + 1
            complex.make_port(port=port, value=qty)
            system.connect((complex, port), species_anchors[species])

    def _build_species(self, system, param_dict):
        species_anchors = {}
        for species, n_reactions in self._species.items():
            this_species = new(
                "Ce", library=LIBRARY, name=species, value=param_dict
            )
            system.add(this_species)

            if n_reactions == 1:
                species_anchors[species] = this_species
            else:
                anchor = new("0", name=species)
                system.add(anchor)
                system.connect(this_species, anchor)
                species_anchors[species] = anchor

            if species in self._flowstats:
                flowstat = new("Sf", value=self._flowstats[species], name=species)
                system.add(flowstat)
                system.connect(flowstat, species_anchors[species])

            if species in self._chemostats:
                chemostat = new("Se", value=0, name=species)
                system.connect(chemostat, species_anchors[species])
                this_species.initial_values['p_0'] = self._chemostats[species]

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
            idx = "{name}".format(name=name) # edit: removed base for named reactions
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

