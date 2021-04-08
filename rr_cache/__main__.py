from rr_cache.rr_cache import (
    rrCache,
)
from rr_cache.Args import (
    add_arguments,
    build_args_parser
)
from logging import (
    Logger
)
from colored import fg, bg, attr
from logging import StreamHandler


def init(
    parser: 'ArgumentParser',
    args: 'Namespace'
) -> Logger:
    from brs_utils import create_logger
    from rr_cache._version import __version__

    if args.log.lower() in ['silent', 'quiet'] or args.silent:
        args.log = 'CRITICAL'

    # Create logger
    logger = create_logger(parser.prog, args.log)

    logger.info(
        '{color}{typo}rr_cache {version}{rst}{color} ({prog}){rst}\n'.format(
            prog = logger.name,
            version = __version__,
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset')
        )
    )
    logger.debug(args)

    return logger


def entry_point():
    parser = build_args_parser(
        prog = 'rr_cache',
        description = 'RetroRules Cache'
    )
    args  = parser.parse_args()

    logger = init(parser, args)

    if args.cache_dir:
        print("rrCache is going to be generated into " + args.cache_dir)
        gen_cache(args.cache_dir, logger)


def gen_cache(outdir, logger):
    rrCache.generate_cache(outdir, logger)


if __name__ == '__main__':
    entry_point()
