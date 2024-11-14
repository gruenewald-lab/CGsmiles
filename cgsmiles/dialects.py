from inspect import signature, Signature, Parameter
from functools import partial

def check_and_cast_types(bound_args, signature):
    for name, value in bound_args.arguments.items():
        param = signature.parameters.get(name)
        # Check if a type annotation is present
        if param and param.annotation != Parameter.empty:
            expected_type = param.annotation

        # Attempt type casting if the value is not of the expected type
        if not isinstance(value, expected_type):
            try:
               bound_args.arguments[name] = expected_type(value)
            except (TypeError, ValueError):
                raise TypeError(f"Argument '{name}' must be of type {expected_type.__name__}")
    return bound_args

def _parse_node(string_iterable,
                dialect_signature,
                annotation_sep_token=';',
                annotation_assign_token='='):
    """
    This base function parsers a CGSmiles node. It must be
    decorated with a signature which defines the dialect.
    The dialect sets expected labels and default values of
    a given node.
    """
    args_found = []
    kwargs_found = {}
    if len(string_iterable) > 0:
        elements = string_iterable.split(annotation_sep_token)
        for entry in elements:
            key_value = entry.split(annotation_assign_token)
            if len(key_value) == 1:
                args_found.append(key_value[0])
            else:
                kwargs_found[key_value[0]] = key_value[1]

    applied_labels = dialect_signature.bind(*args_found,
                                            **kwargs_found)
    applied_labels = check_and_cast_types(applied_labels,
                                          dialect_signature)
    applied_labels.apply_defaults()
    return applied_labels.arguments

def create_dialect(default_attributes):
    """
    Creates a signature of default annotations.
    Note that the order of the entries in the dict
    determines the order of the args accepted.
    """
    parameters = []
    for argname, default_value in default_attributes.items():
        arg_type = type(default_value)
        parameters.append(Parameter(argname,
                                    Parameter.POSITIONAL_OR_KEYWORD,
                                    default=default_value,
                                    annotation=arg_type))
    sig = Signature(parameters)
    return sig

##########################################################
#                   KNOWN DIALECTS                       #
##########################################################
# this one is for global use
# it is the base CGSmiles dialect
GRAPH_BASE = create_dialect({"fragname": "NaN",
                             "c": 0.0,
                             "w": 1.0})
parse_graph_base_node = partial(_parse_node, dialect_signature=GRAPH_BASE)
# this one is an internal fukery until the pysmiles
# base parser is available
# it just strips the kwargs from fragments before
# they go to the respective parser
# in case of cgsmiles fragments it is a bit doing
# double the work
fragment_base = create_dialect({"w": 1.0})
_fragment_node_parser = partial(_parse_node, dialect_signature=fragment_base)
