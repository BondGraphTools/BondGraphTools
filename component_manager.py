import json
import pathlib
import os

LIB_DIR = os.path.join(os.path.dirname(__file__), "components")
libraries = dict()
base_id = "base"

for filename in pathlib.Path(LIB_DIR).glob("**/*.json"):

    with open(filename) as file:
        lib = json.load(file)
        lib_id = lib["id"]
        del lib["id"]
        libraries[lib_id] = lib


def get_library_list():
    return [(l, libraries[l]["name"]) for l in libraries if l != base_id]



def get_components_list(library):

    components = libraries[library]["components"]
    print(components)
    return [(comp_id, components[comp_id]["name"]) for comp_id in components]








