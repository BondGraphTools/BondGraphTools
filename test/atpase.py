import BondGraphTools as bgt
import pathlib
import sys
import os
path = os.path.relpath(path=".")
sys.path.insert(0, os.path.abspath(path))

log_file = "profile.log"


def atpase():

    # Reaction Conductance in fmol/s
    reactions = [
        # reactants, products, rate
        (("P1",), ("Ki", "P2"), 330.5462),
        (("P2",), ("Ki", "P3",), 132850.9145),
        (("P3", "Nai"), ("P4",), 200356.0223),
        (("P4", "Nai"), ("P5",), 2238785.3951),
        (("P5", "Nai"), ("P6", "zF_5"), 10787.9053),
        (("P6",), ("MgADP", "P7"), 15.3533),
        (("P7",), ("P8",), 2.3882),
        (("P8",), ("P9", "Nae", "zF_8"), 2.2855),
        (("P9",), ("P10", "Nae"), 1540.1349),
        (("P10",), ("P11", "Nae"), 259461.6507),
        (("P11", "Ke"), ("P12",), 172042.3334),
        (("P12", "Ke"), ("P13",), 6646440.3909),
        (("P13",), ("P_i", "H", "P14"), 597.4136),
        (("P14", "MgATP"), ("P15",), 70.9823),
        (("P15",), ("P1",), 0.015489)
    ]

    # Chemical Affinities in 1/fmol
    pump_affinities = [101619537.2009,
                       63209.8623,
                       157.2724,
                       14.0748,
                       5.0384,
                       92.6964,
                       4854.5924,
                       15260.9786,
                       13787022.8009,
                       20459.5509,
                       121.4556,
                       3.1436,
                       0.32549,
                       156.3283,
                       1977546.8577]

    species_affinities = {
        "Ki": 0.0012595,  # intracellular potassium
        "Ke": 0.009236,  # extracelluar potassium
        "Nai": 0.00083514,  # intracellular sodium
        "Nae": 0.0061242,  # extracelluar sodium
        "MgATP": 2.3715,
        "MgADP": 7.976e-5,
        "P_i": 0.04565,
        "H": 0.04565
    }

    membrane_capacitance = 153400  # fentoFarad ?
    F = 96485.33289  # faradays constant
    r5_charge = -0.0550 * F
    r8_charge = -0.9450 * F
    Gas_Constant = 8.3144621
    Temperature = 310

    # Initial Conditions
    initial_conditions = {
        "H": 0.0015128,
        "Ke": 27.9828,
        "Ki": 0.038,
        "MgADP": 0.038,
        "MgATP": 380,
        "Nae": 777.3,
        "Nai": 1900,
        "P_i": 0.038,
        "mem": -18408,
        "P1": 4.06727417747628e-10,
        "P2": 1.84495437200145e-09,
        "P3": 0.00183103460985586,
        "P4": 0.0324604777886237,
        "P5": 0.143884405731667,
        "P6": 0.00968432811008776,
        "P7": 0.534480003169437,
        "P8": 0.169644923453282,
        "P9": 5.63940980127252e-07,
        "P10": 7.9739635710452e-05,
        "P11": 0.00282186679243115,
        "P12": 0.0281498805518011,
        "P13": 0.0702590453291214,
        "P14": 0.00626267917464576,
        "P15": 0.000446049400255497
    }
    model = bgt.new(name="ATPase")
    # add all the pumps
    for pump, affinity in enumerate(pump_affinities):
        parameters = {'k': affinity,
                      'R': Gas_Constant,
                      'T': Temperature}
        name = "P{}".format(pump + 1)
        Ce = bgt.new(component="Ce",
                     name=name,
                     library="BioChem",
                     value=parameters)
        junction = bgt.new("0", name=name)

        model.add(Ce, junction)
        model.connect(Ce, junction)

    # add all the species
    for species, affinity in species_affinities.items():
        parameters = {'k': affinity,
                      'R': Gas_Constant,
                      'T': Temperature}
        Ce = bgt.new(component="Ce",
                     name=species,
                     library="BioChem",
                     value=parameters)

        junction = bgt.new("0", name=species)
        model.add(Ce, junction)
        model.connect(Ce, junction)

    # add the electrical components
    c_mem = bgt.new('C', name="mem", value={'C': membrane_capacitance})
    c_junction = bgt.new("0", name="mem")
    model.add(c_mem, c_junction)
    model.connect(c_mem, c_junction)

    zF_5 = bgt.new('TF', name='zF_5', value={'r': r5_charge})
    zF_8 = bgt.new('TF', name='zF_8', value={'r': r8_charge})

    model.add(zF_5, zF_8)
    # port 0 is the primary winding of a transformer and 1 is the secondary
    # winding.
    port = 1
    zF_5_secondary = (zF_5, port)
    zF_8_secondary = (zF_8, port)
    model.connect(zF_5_secondary, c_junction)
    model.connect(zF_8_secondary, c_junction)

    # add all the reactions and wire it up
    for reaction_number, (reactants, products, rate) in enumerate(reactions):
        parameters = {'r': rate, 'R': Gas_Constant, 'T': Temperature}
        name = "R{}".format(reaction_number + 1)

        Re = bgt.new("Re", name=name, library="BioChem", value=parameters)
        model.add(Re)

        if len(reactants) > 1:
            in_port = bgt.new("Y", name=name + "_{in}", library="BioChem")
            model.add(in_port)
            model.connect((Re, 0), (in_port, 0))
            for reactant in reactants:
                junction = model.find(component="0", name=reactant)
                model.connect(junction, in_port)
        else:
            junction = model.find(component="0", name=reactants[0])
            model.connect((Re, 0), junction)

        if len(products) > 1:
            out_port = bgt.new("Y", name=name + "_{out}", library="BioChem")
            model.add(out_port)
            model.connect((Re, 1), (out_port, 0))
            for product in products:
                if product == 'zF_5':
                    junction = (
                        zF_5, 0)  # the primary winding of the charge transporter
                elif product == 'zF_8':
                    junction = (zF_8, 0)
                else:
                    junction = model.find(component="0", name=product)

                model.connect(junction, out_port)
        else:
            junction = model.find(component="0", name=products[0])
            model.connect((Re, 1), junction)

    from math import log
    # Add Chemostats
    for species, affinity in species_affinities.items():
        effort = log(affinity * initial_conditions[species])
        Se = bgt.new('Se', value=effort, name=species)
        model.add(Se)
        junction = model.find(name=species, component='0')
        model.connect(Se, junction)

    # generate equations
    model_equations = model.constitutive_relations
