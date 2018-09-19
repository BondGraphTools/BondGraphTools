"""The file save/load interface and file format data model

This module provides the basic IO functionality such as saving and loading
to file.

todo:
    As the file format matures, we probably want to move to a object
    oriented loader that tracks the (as yet undefined) schema.

"""

import logging
import pathlib
import yaml

from .compound import BondGraph
from .actions import connect, new, expose
from .exceptions import *

logger = logging.getLogger(__name__)

FILE_VERSION = "0.1"

def save(model, filename):
    """
    Saves the model to the specified path

    Args:
        model: The model to be saved
        filename: The file to save to

    """

    model_directory = _build_model_directory(model)

    models = {}
    templates = {}

    for uri, sub_model in model_directory.items():
        models.update({
            uri: _build_model_data(sub_model, templates)
        })
    data = {
        "version":FILE_VERSION,
        "root": model.name,
        "models": models
    }

    with open(filename, 'w') as filestream:
        yaml.dump(data, filestream, default_flow_style=False)


def _build_model_directory(model):

    try:
        directory = {model.uri:model}
        for c in model.components:
            directory.update(_build_model_directory(c))
        return directory
    except AttributeError:
        return {}


def _build_model_data(model, templates):
    components = []
    out = {}
    for c in model.components:
        if isinstance(c, BondGraph):
            components.append(
                f"{c.name} {c.uri}"
            )
        else:
            components.append(
                _build_component_string(c)
            )
    out.update({"components": components})

    netlist = []
    for tail, head in model.bonds:
        netlist.append(f"{tail} {head}")

    if netlist:
        out.update({"netlist": netlist})

    ports = []
    for port in model.ports:
        (c1, e), (c2, f) = model._port_map[port]
        if c1==c2:
            ports.append(f"{c1.name} {port.name}")
        else:
            raise NotImplementedError
    if ports:
        out.update({"ports": ports})

    return out


def _build_component_string(component):

    out_str = f"{component.name} {component.template}"
    logger.debug("Trying to serialise: %s", out_str)
    try:
        for param, value in component.params.items():
            logger.debug("Param: %s, %s", param, value)
            if isinstance(value, (int, float)):
                out_str += f" {param}={value}"
            elif isinstance(value, dict):
                try:
                    v = value["value"]
                    if isinstance(v, (float, int)):
                        out_str += f" {param}={v}"
                except KeyError as ex:
                    logger.debug("Skipping: %s ", str(ex))
                    pass

    except AttributeError:
        pass
    logger.debug("Saving component string: %s", out_str )
    return out_str


def load(file_name, model=None, as_name=None):
    """
    Loads a model from file

    Args:
        file_name (str or Path): The file to load.

    Returns: An instance of `BondGraph`

    """
    if isinstance(file_name, pathlib.Path):
        file_name = str(file_name)

    with open(file_name, 'r') as f:
        data = yaml.load(f)

    version = str(data['version'])
    if version == FILE_VERSION:
        return _builder(data, model, as_name)
    else:
        raise NotImplementedError


def _builder(data, model=None, as_name=None):
    if not model:
        root = "/"
    else:
        root = model
    models = data['models']

    def _build(model_name, template_name):
        model_data = models[template_name]
        netlist = model_data["netlist"]

        logger.debug("%s: Trying to build", model_name)

        model = new(name=model_name)

        for comp_string in model_data["components"]:
            logger.debug("%s: building", comp_string)

            try:
                comp = _base_component_from_str(comp_string)
            except (ValueError, KeyError):
                name, sub_model = comp_string.split(" ")
                comp = _build(name, sub_model)

            model.add(comp)
            logger.debug("%s components complete", model_name)

        _wire(model, netlist)

        try:
            IO_ports = model_data["ports"]
            _expose(model, IO_ports)
        except KeyError:
            logger.debug("No ports on model ")

        return model

    def _wire(model, netlist):
        logger.debug("%s: trying to wire", model.name)

        def get_port(port_string):

            tokens = iter(port_string.split('.'))

            c = next(tokens)
            try:
                comp, = (comp for comp in model.components if comp.name == c)
            except ValueError:
                raise InvalidComponentException(
                    f"Could not find component {c} in model {model.name}")
            try:
                t2 = next(tokens)
            except StopIteration:
                return comp
            try:
                t3 = next(tokens)
            except StopIteration:

                logger.debug("Tyring to get port %s, %s", str(comp), str(t2))
                try:
                    t2 = int(t2)
                except ValueError:
                    pass
                port = comp.get_port(t2)
                logger.debug("Got %s", str(port))
                return port

            else:
                raise NotImplementedError

        for bond_string in netlist:
            logger.debug("%s: bond %s", model.name, bond_string)
            tail_str, head_str = bond_string.split()
            tail = get_port(tail_str)
            head = get_port(head_str)
            connect(tail, head)

    def _expose(model, IO_ports):
        for port_string in IO_ports:
            component_name, port_label = port_string.split(" ")
            comp, = {
                c for c in model.components if c.name == component_name
            }
            expose(comp, port_label)

    def _parse_build_args(in_args):

        if not in_args:
            return [] , {}

        arg = in_args[0]
        args, kwargs = _parse_build_args(in_args[1:])

        try:
            k, v = arg.split("=")
        except ValueError:
            k = None
            v = arg

        if v.isnumeric():
            try:
                v = int(v)
            except ValueError:
                v= float(v)
        if not k:
            args.append(v)
        else:
            kwargs.update({k:v})

        return args, kwargs

    def _base_component_from_str(string):
        label, tempate, *build_args = string.split()
        args, kwargs = _parse_build_args(build_args)
        library, component = tempate.split("/")

        comp = new(name=label,
                   library=library,
                   component=component,
                   value=args)

        for k, v in kwargs.items():
            comp.set_param(k, v)
        return comp

    def _validate(model_dict):
        # We should use this function to make sure we don't have an infinite loop
        # in the build cycle.
        pass

    _validate(models)

    out = _build(root, root)

    if as_name:
        out.name = as_name
    elif not model:
        out.name = data["root"]

    return out



