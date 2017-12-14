import json
import pathlib
import os

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
    return [(l, __libraries[l]["desc"]) for l in __libraries if l != base_id]


def get_components_list(library):

    components = __libraries[library]["components"]

    return [(comp_id, components[comp_id]["desc"]) for comp_id in components]


def get_component_data(library, component):

    return __libraries[library]["components"][component]


def find(component, restrict_to=None, find_all=False):

    results = []

    if restrict_to:
        keys = {k for k in __libraries if k in restrict_to}
    else:
        keys = __libraries

    for lib_id in keys:
        for comp_id in __libraries[lib_id]["components"]:
            if comp_id == component:
                if find_all:
                    results.append((lib_id, comp_id))

                else:
                    return lib_id, comp_id
    if find_all:
        return results
    else:
        raise NotImplementedError("Component not found", component)




