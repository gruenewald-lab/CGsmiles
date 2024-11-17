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

def _parse_dialect_string(string_iterable,
                          dialect_signature,
                          annotation_sep_token=';',
                          annotation_assign_token='='):
    """
    This base function parsers a string that describes key value pairs
    in having a pattern of:

    key<annotation_assign_token>value<annotation_sep_token>key ...

    Default values, non-keyword agruments and types are defined using the
    dialect signature object. If args are defined the key and assignment
    token may be omitted.

    Neither the `annotation_sep_token` nor the `annotation_assign_token`
    can be part of key or value. A SyntaxError is raised in this case.

    Parameters
    ----------
    string_iterable: iter
        the string or iter object that contains the string
    dialect_signature: cls.inspec.Signature
        a signature defineing args, kwargs, default values
        and types
    annotation_sep_token: str
        character used to seperate key value pairs
    annotation_assign_token: str
        character used to assign a key from a value

    Returns
    -------
    dict
       dict of key value paris

    Raises
    ------
    SyntaxError
        an error is raised if the signature does not match or
        too many annotation_assign_token are given
    """
    args_found = []
    kwargs_found = {}
    if len(string_iterable) > 0:
        elements = string_iterable.split(annotation_sep_token)
        for entry in elements:
            if entry.count('=') > 1:
                # this takes care of too many '=' chacaters
                msg = (f"Your annotation {entry} contains too many "
                       f"{annotation_assign_token} charachters. Only"
                        "chacracter per key value pair is allowed")
                raise SyntaxError(msg)
            key_value = entry.split(annotation_assign_token)

            if len(key_value) == 1:
                args_found.append(key_value[0])
            else:
                kwargs_found[key_value[0]] = key_value[1]

    try:
        applied_labels = dialect_signature.bind(*args_found,
                                                **kwargs_found)
    except TypeError as emsg:
        print(emsg)
        msg = ("You have too many positional arguments or "
               f"{annotation_sep_token} as part of key value "
                "pairs which is not allowed.")
        raise SyntaxError(msg)

    applied_labels = check_and_cast_types(applied_labels,
                                          dialect_signature)
    applied_labels.apply_defaults()
    # if there are kwargs we need to put them into
    # output dict
    out_args = {}
    out_args.update(applied_labels.arguments['kwargs'])
    del applied_labels.arguments['kwargs']
    out_args.update(applied_labels.arguments)
    return out_args

def create_dialect(default_attributes, accept_kwargs=True):
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
    if accept_kwargs:
        parameters.append(Parameter('kwargs',
                                    kind=Parameter.VAR_KEYWORD))
    sig = Signature(parameters)
    return sig

##########################################################
#                   KNOWN DIALECTS                       #
##########################################################
# this one is for global use
# it is the base CGSmiles dialect
CGSMILES_DEFAULT_DIALECT = create_dialect({"fragname": "NaN",
                                           "c": 0.0,
                                           "w": 1.0})
parse_graph_base_node = partial(_parse_dialect_string, dialect_signature=CGSMILES_DEFAULT_DIALECT)
# this one is an internal fukery until the pysmiles
# base parser is available
# it just strips the kwargs from fragments before
# they go to the respective parser
# in case of cgsmiles fragments it is a bit doing
# double the work
fragment_base = create_dialect({"w": 1.0})
_fragment_node_parser = partial(_parse_dialect_string, dialect_signature=fragment_base)
