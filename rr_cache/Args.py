from argparse  import ArgumentParser
from typing import Callable
from os import path as os_path

from brs_utils import add_logger_args


HERE = os_path.dirname(os_path.abspath( __file__ ))
DATA_PATH = os_path.join(HERE, 'data')
# Default values for the arguments
DEFAULTS = {
    'data_type': 'mnx3.1',
    'attrs': [],
    'interactive': False,
    'do_not_dwnl_cache': False,
}


def build_args_parser(
    prog: str,
    description: str = '',
    epilog: str = '',
    m_add_args: Callable = None,
) -> ArgumentParser:

    parser = ArgumentParser(
        prog = prog,
        description = description,
        epilog = epilog
    )

    # Build Parser with rptools common arguments
    parser = add_arguments(parser)

    # Add module specific arguments
    if m_add_args is not None:
        parser = m_add_args(parser)

    return parser


def add_arguments(parser: ArgumentParser) -> ArgumentParser:
    # Add arguments related to the logger
    parser = add_logger_args(parser)

    parser.add_argument(
        '--data-type',
        default=DEFAULTS['data_type'],
        type=str,
        help='Type of data to use (e.g. mnx3.1, mnx4.0...). Determines which configuration files and folders to use both the cache and the input cache (default: %(default)s).'
    )
    parser.add_argument(
        '--gen-cache',
        default=None,
        action='store_true',
        help='generate the cache and exits'
    )
    parser.add_argument(
        '--reaction-rules',
        default=None,
        type=str,
        nargs='*',
        help='print out specified reaction rule(s), all if nothing specified'
    )
    parser.add_argument(
        '--reactions',
        default=None,
        type=str,
        nargs='*',
        help='print out specified reaction(s), all if nothing specified'
    )
    parser.add_argument(
        '--compounds',
        default=None,
        type=str,
        nargs='*',
        help='print out specified compound(s), all if nothing specified'
    )
    parser.add_argument(
        '--interactive',
        default=DEFAULTS['interactive'],
        action='store_true',
        dest='interactive',
        help='Ask user for confirmation when loading cache'
    )
    parser.add_argument(
        '--do-not-dwnl-cache',
        default=DEFAULTS['do_not_dwnl_cache'],
        action='store_true',
        help='Do not download the cache from the remote repository'
    )

    return parser
