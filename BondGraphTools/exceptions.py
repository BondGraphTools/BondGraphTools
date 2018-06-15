import logging
import logging.config


class InvalidPortException(Exception):
    pass


class InvalidComponentException(Exception):
    pass


class ModelParsingError(Exception):
    pass


class ModelException(Exception):
    pass


class SymbolicException(Exception):
    pass
