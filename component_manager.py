"""
Component Library Manager
=========================

This module takes care of the loading and management of component libraries.
This library will automatically load all factory libraries upon import.
Additionally libraries can be added via :func:`load_library`.

Component Libraries are expected to be in json format.
The structure must be::
    {
        "id": The unique library id (str),
        "desc": The description of this library
        "components": {
            "component_id": component dictionary,
            ...
        }
    }

Each Component dictionary must be of the form::
    {
        "desc": Human friendly description (str)
        "cls": ("OnePort", "TwoPort", "IOPort", "ManyPort")
        "string": (optional) Overrides display string
        "ports": See Note,
        "local_params": (optional) See Note
        "global_params": (optional) See Note
        "constitutive_relations": See Note
    }
"""

import json
import pathlib
import os
import logging

logger = logging.getLogger(__name__)


LIB_DIR = os.path.join(os.path.dirname(__file__), "components")
__libraries = dict()
base_id = "base"


def load_library(filename):
    """
    Loads the library specified by the filename and makes it available for use.

    Args:
        filename: The (absolute) filename of the library to be loaded.

    Returns:
        bool: True if the library was successfully loaded.
    """
    with open(filename) as file:
        try:
            lib = json.load(file)
            lib_id = lib["id"]
            del lib["id"]
            if lib_id in __libraries:
                raise KeyError("Invalid Library ID")
            elif set(lib.keys()) != {"desc", "components"}:
                raise KeyError("Invalid Library")
            __libraries[lib_id] = lib
        except json.JSONDecodeError as ex:
            logger.warning("Error loading library %s; %s",
                           filename,
                           ex.args)
            return False
        except KeyError as ex:
            logger.warning("Error loading library %s: %s",
                           filename, ex.args[0])
            return False
    return True


def get_library_list():
    """
    Fetches a list of the libraries available for use.

    Returns:
        list of (library id, description) tuples
    """

    return [(l, __libraries[l]["desc"]) for l in __libraries if l != base_id]


def get_components_list(library):
    """
    Fetches a list of components available in the given library.

    Args:
        library: The library id of the library to query

    Returns:
        list of (component id, description) tuples
    """
    components = __libraries[library]["components"]

    return [(comp_id, components[comp_id]["desc"]) for comp_id in components]


def get_component_data(library, component):
    """
    Fetches the component data for the specified component
    Args:
        library: The id of the library to which the component belongs
        component: The id of the specific component

    Returns:
        dict - the component dictionary

    """
    return __libraries[library]["components"][component]


def find(component, restrict_to=None, find_all=False, ensure_unique=False):
    """
    Finds the specified component.

    Args:
        component: The component id to find.
        restrict_to: `list` or `set` of library id's to be that the search
         should be restricted to.
        find_all: `False` if the function should return only the first instance
         of the component, `True` if the function should return all such
         instances
         ensure_unique: If true, this assumes that the component id must be
          unique across libraries, and hence will raise an exception if this is
          assumption is violated.

    Returns: the library id, or a list of library_id in which this component
     can be found.

    Raises:
        NotImplementedError - if the component is not found.
        ValueError - if the component is assume to be unique, but is not
    """
    results = []

    if restrict_to:
        keys = {k for k in __libraries if k in restrict_to}
    else:
        keys = __libraries
    unique_id = None
    for lib_id in keys:
        for comp_id in __libraries[lib_id]["components"]:
            if comp_id == component:
                if find_all:
                    results.append(lib_id)
                elif not ensure_unique:
                    return lib_id
                elif ensure_unique and not unique_id:
                    unique_id = lib_id
                elif ensure_unique and unique_id:
                    raise ValueError(
                        "Could not find component: Component not unique",
                        component)
    if find_all:
        return results
    elif unique_id:
        return unique_id
    else:
        raise NotImplementedError("Component not found", component)


## On Import actions
for filename in pathlib.Path(LIB_DIR).glob("**/*.json"):
    load_library(filename)