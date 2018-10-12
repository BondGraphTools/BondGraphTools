"""Exceptions and errors for BondGraphTools"""


class InvalidPortException(Exception):
    """Exception for trying to access a port that is in use, or does not
    exist """
    pass


class InvalidComponentException(Exception):
    """Exception for when trying to use a model that can't be found, or is of
    the wrong type"""


class ModelParsingError(Exception):
    """Exception for problems generating symbolic equations from string"""


class ModelException(Exception):
    """Exception for inconsistent or invalid models when running simulations"""


class SymbolicException(Exception):
    """Exception for when there are issues in model reduction or symbolic
    manipulation"""


class SolverException(Exception):
    """Exception for issues running numerical solving."""
