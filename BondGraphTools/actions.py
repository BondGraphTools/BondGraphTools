"""This module provides functions for the actions one wishes to perform on
models bond graph models such as creating new models and components,
connecting parts together, and setting parameters.
"""

import copy
import logging

from BondGraphTools.component_manager import get_component, base_id
from BondGraphTools.exceptions import (
    InvalidPortException, InvalidComponentException
)
from BondGraphTools.base import BondGraphBase, Bond, Port
from BondGraphTools.port_managers import PortTemplate

from .atomic import EqualFlow
# from .port_hamiltonian import PortHamiltonian

logger = logging.getLogger(__name__)

__all__ = [
    "new",
    "add",
    "swap",
    "connect",
    "disconnect",
    "expose",
    "remove",
    "set_param"
]


def disconnect(target, other):
    """
    Disconnects the flow of energy between the two components or ports.
    If there is no connection, this method does nothing.

    Args:
        target (Port, BondGraphBase):
        other (Port, BondGraphBase):

    Raises:
        InvalidComponentException

    See Also:
        :func:`connect`
    """
    if isinstance(target, BondGraphBase):
        model = target.parent
    elif isinstance(target, Port):
        model =  target.component.parent
    else:
        model = target[0].parent

    if isinstance(other, BondGraphBase):
        model_prime = other.parent
    elif isinstance(other, Port):
        model_prime = other.component.parent
    else:
        model_prime = other[0].parent

    if not model or not model_prime or (model is not model_prime):
        raise InvalidComponentException(f"Could not find components")

    def _filter(item):
        # assume item is a port:
        if isinstance(item, Port):
            return {bond for bond in model.bonds if item in bond}
        try:
            _, _ = item
            return {bond for bond in model.bonds
                    if item in (bond.head, bond.tail)}
        except TypeError as ex:
            return {bond for bond in model.bonds
                    if item is bond.head.component or
                    item is bond.tail.component}

    targets = _filter(target) & _filter(other)

    for target_bond in targets:
        model._bonds.remove(target_bond)


def connect(source, destination):
    """Connects two components or ports.

    Defines a power bond between the source and destination ports such that
    the bond tail is at the source, and the bond head is at the destination.
    We assume that either the source and/or destination is or has a free
    port.

    Args:
        source (Port or BondGraphBase):      The tail of the power bond
        destination (Port or BondGraphBase): The head of the power bond

    Raises:
        InvalidPortException, InvalidComponentException

    See Also:
        :func:`disconnect`
    """


    tail = _find_or_make_port(source, is_tail=True)
    head = _find_or_make_port(destination)

    model = tail.component.parent
    if not model or head.component.parent is not model:
        raise InvalidComponentException
    bond = Bond(tail, head)
    model._bonds.add(bond)


def _unpack_port_arg(port):
    # returns (component, idx)
    if isinstance(port, (Port, tuple)):
        c, idx = port
        return c, idx

    if isinstance(port, BondGraphBase):
        return port, None

    if isinstance(port, PortTemplate):
        assert port.parent
        return port.parent, port

    raise InvalidPortException("Could not unpack port %s", port)


def _find_or_make_port(arg, is_tail=False):
    component, port = _unpack_port_arg(arg)

    if isinstance(component, EqualFlow) and port is None:
        if is_tail:
            port = component.inverting
        else:
            port = component.non_inverting
    try:
        p = component.get_port(port)
        if p is None:
            raise InvalidPortException
        return p
    except InvalidPortException as ex:
        try:
            return component.new_port(port)
        except AttributeError:
            raise ex


def swap(old_component, new_component):
    """
    Replaces the old component with a new component.
    Components must be of compatible classes; 1 one port cannot replace an
    n-port, for example.
    The old component will be completely removed from the system model.

    Args:
        old_component: The component to be replaced. Must already be in the
         model.
        new_component: The substitute component which must not be in the
         model

    Raises:
        InvalidPortException, InvalidComponentException
    """

    # TODO: More validation required
    #
    def is_swap_valid(old_comp, new_comp):
        if not isinstance(new_comp, BondGraphBase):
            return False

        # Dirty Hack because 'in' for BondSet double counts
        num_bonds = len({b for b in old_comp.parent.bonds if old_comp is
                         b.head.component or old_comp is b.tail.component})

        # if isinstance(new_comp, PortManager) and len(new_comp.ports) < num_bonds:
        #     return False

        return True

    model = old_component.parent
    if not model:
        return

    if new_component in model.components:
        raise InvalidComponentException(
            "Component is already in the model"
        )

    elif not is_swap_valid(old_component, new_component):
        raise InvalidComponentException("Cannot swap components")

    model.add(new_component)

    swaps = []
    for bond in model.bonds:
        try:
            if bond.tail.component is old_component:
                tail = new_component.get_port()
                head = bond.head
            elif bond.head.component is old_component:
                tail = bond.tail
                head = new_component.get_port()
            else:
                continue
        except InvalidPortException:
            raise InvalidComponentException(
                "Cannot swap components: Incompatible ports"
            )

        swaps.append((bond, Bond(tail, head)))

    for old_bond, new_bond in swaps:
        disconnect(old_bond.tail, old_bond.head)
        connect(new_bond.tail, new_bond.head)

    model.remove(old_component)


def new(component=None, name=None, library=base_id, value=None, **kwargs):
    """
    Creates a new Bond Graph from a library component.

    Args:
        component(str or obj): The type of component to create.
         If a string is specified, the the component will be created from the
         appropriate libaray. If an existing bond graph is given, the bond
         graph will be cloned.
        name (str): The name for the new component
        library (str): The library from which to find this component (if
        component is specified by string).
        value:

    Returns: instance of :obj:`BondGraph`

    Raises: NotImplementedError
    """

    if not component:
        cls = _find_subclass("BondGraph", BondGraphBase)
        return cls(name=name)
    elif isinstance(component, str):
        build_args = get_component(component, library)
        args = ()
        if name:
            build_args.update({"name": name})
        if value or isinstance(value, (int, float, complex)):
            args, build_args = _update_build_params(args,
                                                    build_args,
                                                    value,
                                                    **kwargs)
        cls = _find_subclass(
            build_args["class"], BondGraphBase
        )
        del build_args["class"]

        comp = cls(*args, **build_args)
        comp.__component__ = component
        comp.__library__ = library
        return comp

    elif isinstance(component, BondGraphBase):
        obj = copy.copy(component)
        if name:
            obj.name = name
        if value:
            _update_build_params(obj.__dict__, value)

        return obj

    else:
        raise NotImplementedError(
            "New not implemented for object {}", component
        )


def _update_build_params(args, build_args, value, **kwargs):

    if isinstance(value, (list, tuple)):
        assignments = zip(build_args["params"].keys(), value)

        for param, v in assignments:
            try:
                build_args["params"][param]["value"] = v
            except TypeError:
                build_args["params"][param] = v

    elif isinstance(value, dict):
        for param, v in value.items():
            try:
                if isinstance(build_args["params"][param], dict):
                    build_args["params"][param]["value"] = v
                else:
                    build_args["params"][param] = v
            except KeyError:
                build_args[param] = v

    elif isinstance(value, str):
        args = (*args, value)
    else:
        # Todo: fix me! Dirty Hacks to make PH load

        p = next(iter(build_args["params"]))
        build_args["params"][p] = value

    return args, build_args


def _find_subclass(name, base_class):

    for c in base_class.__subclasses__():
        if c.__name__ == name:
            return c
        else:
            sc = _find_subclass(name, c)
            if sc:
                return sc


def expose(component, label=None):
    """
    Exposes the component as port on the parent.

    If the target component is not a SS component, it is replaced with a new
    SS component.
    A new external port is added to the parent model, and connected to the
    SS component.

    Args:
        component: The component to expose.
        label: The label to assign to the external port

    Raises: InvalidComponentException
    """
    model = component.parent
    if not model:
        raise InvalidComponentException(
            f"Component {component} is not inside anything"
        )
    # fix me with metamodeles or something trickier
    if component.__component__ is not "SS":
        ss = new("SS",  name=component.name)
        try:
            swap(component, ss)
        except InvalidComponentException as ex:
            raise InvalidComponentException(f"Cannot expose {component}", ex)
    else:
        ss = component

    effort_port = (ss, 'e')
    flow_port = (ss, 'f')

    if not label:
        label = str(len(model.ports))

    model.map_port(label, (effort_port, flow_port))


def add(model, *args):
    """
    Add the specified component(s) to the model
    """
    model.add(*args)


def remove(model, component):
    """Removes the specified components from the Bond Graph model.
    """
    model.remove(component)


def set_param(component, param, value):
    """
    Sets the specified parameter to a particular value.

    Args:
        component (`BondGraphBase`): The particular component.
        param: The parameter to set
        value: The value to assign it to, may be None
    """
    component.set_param(param, value)
