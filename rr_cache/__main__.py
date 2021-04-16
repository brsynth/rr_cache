from rr_cache.rr_cache import (
    rrCache,
)
from rr_cache.Args import (
    add_arguments,
    build_args_parser
)
from logging import (
    Logger,
    getLogger
)
from colored import fg, bg, attr
from logging import StreamHandler
from argparse import (
    ArgumentParser,
    Namespace
)
from json import dumps
from typing import (
    List,
    Dict,
    Tuple
)


def init(
    parser: 'ArgumentParser',
    args: 'Namespace'
) -> Logger:
    from brs_utils import create_logger
    from rr_cache._version import __version__

    if args.log.lower() in ['silent', 'quiet'] or args.silent:
        args.log = 'CRITICAL'

    if args.log.lower() in ['silent', 'quiet', 'def_info'] or args.silent:
        disable_rdkit_logging()
        # # Disable RDKIT logging
        # from rdkit import RDLogger
        # RDLogger.DisableLog('rdApp.*')


    # Create logger
    logger = create_logger(parser.prog, args.log)

    logger.info(
        '{color}{typo}rr_cache {version}{rst}{color}{rst}\n'.format(
            prog = logger.name,
            version = __version__,
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset')
        )
    )
    logger.debug(args)

    return logger


def disable_rdkit_logging():
    """
    Disables RDKit whiny logging.
    """
    import rdkit.rdBase as rkrb
    import rdkit.RDLogger as rkl
    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog('rdApp.error')


def entry_point():
    parser = build_args_parser(
        prog = 'rr_cache',
        description = 'RetroRules Cache'
    )
    args  = parser.parse_args()

    logger = init(parser, args)

    cache = rrCache(
        db=args.db,
        attrs=None,
        logger=logger
    )

    if args.cache_dir:
        # print("rrCache is going to be generated into " + args.cache_dir)
        gen_cache(
            args.cache_dir,
            logger
        )
    elif args.reaction_rules is not None:
        print_rr(
            cache,
            args.reaction_rules,
            logger
        )
    else:
        cache.load(args.attrs)


def gen_cache(outdir, logger):
    rrCache.generate_cache(outdir, logger)


def print_rr(
    cache: 'rrCache',
    reaction_rules: List,
    logger: Logger=getLogger(__file__)
) -> None:
    cache.load(['rr_reactions'])
    if reaction_rules == []:
        print(
            dumps(
                cache.get('rr_reactions'),
                indent=4
            )
        )
    else:
        for rr_id in reaction_rules:
            try:
                print(
                    rr_id+':',
                    dumps(
                        cache.get('rr_reactions')[rr_id],
                        indent=4
                    )
                )
            except KeyError:
                logger.error('Reaction rule ID not found: '+rr_id)


if __name__ == '__main__':
    entry_point()
