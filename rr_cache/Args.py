from argparse  import ArgumentParser
from rr_cache._version import __version__
from typing import(
    Callable,
)


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
    parser.add_argument(
        '--gen-cache',
        default=None,
        type=str,
        dest='cache_dir',
        help='generate the cache and exits'
    )
    parser.add_argument(
        '--reaction-rules',
        default=None,
        type=str,
        nargs='*',
        help='generate the cache and exits'
    )
    parser.add_argument(
        '--db',
        type=str,
        default='file',
        help='Name of the server which stores the cache. If \'file\' is given (default), then cache is loaded from file into memory.'
    )
    parser.add_argument(
        '--attrs',
        type=str,
        choices=[
            'deprecatedCID_cid', 'deprecatedRID_rid', 'deprecatedCompID_compid',
            'cid_strc', 'cid_name', 'cid_xref', 'chebi_cid', 'inchikey_cid',
            'rr_reactions', 'rr_full_reactions',
            'comp_xref',
            'all'
        ],
        default=[],
        nargs='+',
        help='Name(s) of attributes to load, all if not given (default).'
    )
    parser.add_argument(
        '--log', '-l',
        metavar='ARG',
        type=str,
        choices=[
            'debug', 'info', 'warning', 'error', 'critical', 'silent', 'quiet',
            'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', 'SILENT', 'QUIET'
        ],
        default='def_info',
        help='Adds a console logger for the specified level (default: error)'
    )
    parser.add_argument(
        '--silent', '-s',
        action='store_true',
        default=False,
        help='run %(prog)s silently'
    )
    parser.add_argument(
        '--version', '-v',
        action='version',
        version='%(prog)s {}'.format(__version__),
        help='show the version number and exit'
    )
    return parser
