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

from .base import Port, Bond
from .actions import connect, new
from .exceptions import *
logger = logging.getLogger(__name__)

def save(model, filename):
    """
    Saves the model to the specified path

    Args:
        model: The model to be saved
        filename: The file to save to

    """
    raise NotImplementedError

def load(file_name):
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
    if version == "0.1":
        return _builder(data)
    else:
        raise NotImplementedError

def _builder(data):

    def _build(model_name, template_name, models):
        model_data = models[template_name]
        netlist = model_data["netlist"]
        model = new(name=model_name)

        for comp_string in model_data["components"]:
            label, tempate, *build_args = comp_string.split()
            args, kwargs = _parse_build_args(build_args)
            library, component = tempate.split("/")
            try:
                name = kwargs.pop("name")
            except KeyError:
                name = None

            comp = new(name=label,
                    library=library,
                    component=component,
                    value=args)

            for k,v in kwargs.items():
                comp.set_param(k,v)
            model.add(comp)

        _wire(model, netlist)

        return model


    def _wire(model, netlist):
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
                print(t2)
                port = comp.get_port(t2)
                assert port
                return port

            else:
                raise NotImplementedError

        for bond_string in netlist:
            tail_str, head_str = bond_string.split()
            tail = get_port(tail_str)
            head = get_port(head_str)
            connect(tail, head)

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

    def _build_component(name, template, *args, **kwargs):
        pass


    def _validate(model_dict):
        # We should use this function to make sure we don't have an infinite loop
        # in the build cycle.
        pass

    root = data["root"]
    models = data['models']
    _validate(models)

    return _build(root,root, models)



