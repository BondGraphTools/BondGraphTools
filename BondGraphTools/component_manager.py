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
        "description": The description of this library
        "components": {
            "component_id": component dictionary,
            ...
        }
    }

Each Component dictionary must be of the form::
    {

    }
"""
import copy
import json
import pathlib
import os
import logging

logger = logging.getLogger(__name__)

__dir, _ = os.path.split(__file__)
__LIB_DIR = os.path.join(__dir, "components")
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
            elif set(lib.keys()) != {"description", "components"}:
                raise KeyError("Invalid Library")
            __libraries[lib_id] = lib
        except json.JSONDecodeError as ex:
            logger.critical("Error loading library %s; %s",
                            filename,
                            ex.args)
            return False
        except KeyError as ex:
            logger.critical("Error loading library %s: %s",
                            filename, ex.args[0])
            return False

    return True


def get_library_list():
    """
    Fetches a list of the libraries available for use.

    Returns:
        list of (library id, description) tuples
    """

    return [(l, __libraries[l]["description"]) for l in __libraries]


def get_components_list(library):
    """
    Fetches a list of components available in the given library.

    Args:
        library: The library id of the library to query

    Returns:
        list of (component id, description) tuples
    """
    components = __libraries[library]["components"]

    return [(comp_id, components[comp_id]["description"]) for comp_id in components]


def get_component(component, library=base_id):
    """
    Fetches the component data for the specified component
    Args:
        component: The id of the specific component
        library: The id of the library to which the component belongs

    Returns:
        dict - the component dictionary

    """
    return copy.deepcopy(__libraries[library]["components"][component])


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

    results = {l for l in keys if component in __libraries[l]["components"]}

    if find_all:
        return results
    elif (ensure_unique and len(results) == 1) or not ensure_unique:
        return results[0]
    elif ensure_unique:
        raise KeyError("Component is not unique")
    else:
        raise NotImplementedError("Component not found", component)


# On Import actions
for lib_file in pathlib.Path(__LIB_DIR).glob("**/*.json"):
    load_library(lib_file)
