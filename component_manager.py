import json
import pathlib
import os

import model as mo

LIB_DIR = os.path.join(os.path.dirname(__file__), "components")
__libraries = dict()
base_id = "base"

for filename in pathlib.Path(LIB_DIR).glob("**/*.json"):

    with open(filename) as file:
        lib = json.load(file)
        lib_id = lib["id"]
        del lib["id"]
        __libraries[lib_id] = lib


def get_library_list():
    return [(l, __libraries[l]["name"]) for l in __libraries if l != base_id]


def get_components_list(library):

    components = __libraries[library]["components"]

    return [(comp_id, components[comp_id]["name"]) for comp_id in components]


def build_node(library, component):
    kwargs = __libraries[library]["components"][component]

    node = mo.AtomicNode(**kwargs, node_type=component)

    return node






