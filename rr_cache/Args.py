from argparse  import ArgumentParser
from typing import Callable
from os import path as os_path

from brs_utils import add_logger_args


DEFAULTS = {
    'mnx_version': '4.4',
    'cache_dir': '',
    'attrs': []
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
        '--mnx-version',
        default=DEFAULTS['mnx_version'],
        type=str,
        help='Version of MetanetX to use'
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
        '--cache-dir',
        default=DEFAULTS['cache_dir'],
        type=str,
        help='Path to the cache to generate or read from'
    )
    parser.add_argument(
        '--attrs',
        type=str,
        choices=[
            'deprecatedCID_cid', 'deprecatedRID_rid', 'deprecatedCompID_compid',
            'cid_strc', 'cid_name', 'cid_xref', 'chebi_cid', 'inchikey_cid',
            'rr_reactions', 'template_reactions',
            'comp_xref',
            'all'
        ],
        default=DEFAULTS['attrs'],
        nargs='+',
        help='Name(s) of attributes to load, all if not given (default).'
    )

    return parser
