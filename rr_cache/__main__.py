from rr_cache.rr_cache import (
    rrCache,
)
from rr_cache.Args import (
    build_args_parser
)
from logging import (
    Logger,
    getLogger
)
from colored import (
    fg,
    attr
)
from argparse import (
    ArgumentParser,
    Namespace
)
from json import dumps
from typing import (
    List,
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
    logger = create_logger(parser.prog, args.log, 'rr_cache.log')

    logger.info(
        '{color}{typo}rr_cache {version}{rst}{color}{rst}\n'.format(
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
        data_type=args.data_type,
        interactive=args.interactive,
        do_not_dwnl_cache=args.do_not_dwnl_cache,
        logger=logger
    )

    # try:
    if args.gen_cache:
        cache.generate_cache(
            interactive=args.interactive,
            logger=logger
        )
    elif args.reaction_rules is not None:
        print_attr(
            cache,
            'rr_reactions',
            args.reaction_rules,
            args.do_not_dwnl_cache,
            logger
        )
    elif args.reactions is not None:
        print_attr(
            cache,
            'template_reactions',
            args.reactions,
            args.do_not_dwnl_cache,
            logger
        )
    elif args.compounds is not None:
        print_attr(
            cache,
            'cid_strc',
            args.compounds,
            args.do_not_dwnl_cache,
            logger
        )
    else:
        cache.load(interactive=args.interactive, do_not_dwnl_cache=args.do_not_dwnl_cache)
    # except Exception as e:
    #     logger.debug(f"Exception type: {type(e).__name__}")
    #     logger.error(
    #         '\n*** An error occurred:\n{error}'.format(error=str(e))
    #     )
    #     logger.error('\nExiting...\n')
    #     exit(1)

def print_attr(
    cache: 'rrCache',
    attr: str,
    attr_lst: List,
    do_not_dwnl_cache: bool,
    logger: Logger = getLogger(__file__)
) -> None:
    cache.load(attrs=[attr], do_not_dwnl_cache=do_not_dwnl_cache)
    if attr_lst == []:
        print(
            dumps(
                cache.get(attr),
                indent=4
            )
        )
    else:
        for id in attr_lst:
            try:
                print(
                    id+':',
                    dumps(
                        cache.get(attr)[id],
                        indent=4
                    )
                )
            except KeyError:
                logger.error(
                    'ID not found in rrCache(\'{attr}\'): {id}'.format(
                        attr=attr,
                        id=id
                    )
                )


if __name__ == '__main__':
    entry_point()
