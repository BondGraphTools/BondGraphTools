import logging
import pathlib
import yaml

from .base import new

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

    Returns:

    """
    if isinstance(file_name, pathlib.Path):
        file_name = str(file_name)

    with open(file_name, 'r') as f:
        data = yaml.load(f)

    root = data["root"]
    models = data['models']
    _validate(models)
    model = _build(root,root, models)
    return model

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

        comp = new(name=name,
                library=library,
                component=component,
                value=args)

        for k,v in kwargs.items():
            comp.set_param(k,v)

        model.add(comp, label=label)

    for bond_string in netlist:
        c1_str, c2_str = bond_string.split()
        try:
            c1, p1 = c1_str.split('.')
            p1 = int(p1)
        except ValueError:
            c1 = c1_str
            p1 = None
        try:
            c2, p2 = c2_str.split('.')
            p2 = int(p2)
        except ValueError:
            c2 = c2_str
            p2 = None
        model.connect((c1,p1), (c2, p2))
    return model


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




