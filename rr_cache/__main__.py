from rr_cache.rr_cache import (
    rrCache,
)
from rr_cache.Args import add_arguments
from brs_utils import init as init_logger, build_args_parser
from logging import Logger, getLogger
from colored import fg, attr
from json import dumps
from typing import (
    List,
)
from .Args import CONFIG_PATH
from ._version import __version__


def disable_rdkit_logging():
    """
    Disables RDKit whiny logging.
    """
    import rdkit.rdBase as rkrb
    import rdkit.RDLogger as rkl

    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog("rdApp.error")


def entry_point():
    parser = build_args_parser(
        prog="rr_cache",
        version=__version__,
        description="RetroRules Cache",
        m_add_args=add_arguments,
    )
    args = parser.parse_args()
    if args.log.lower() in ["silent", "quiet", "def_info"] or args.silent:
        disable_rdkit_logging()

    logger = init_logger(parser, args, __version__)

    if args.list_chemical_spaces:
        # list config_*.json files in CONFIG_PATH
        import os

        files = os.listdir(CONFIG_PATH)
        cspaces = []
        for file in files:
            if file.startswith("config_") and file.endswith(".json"):
                cspaces.append(file[len("config_") : -len(".json")])
        print(
            "{color}{typo}Available chemical spaces:{rst}{color}{rst}\n".format(
                color=fg("white"), typo=attr("bold"), rst=attr("reset")
            )
        )
        for dt in cspaces:
            print(f"- {dt}")
        print()
        exit(0)

    cache = rrCache(
        cspace=args.cspace,
        databases=args.db,
        interactive=args.interactive,
        do_not_dwnl_cache=args.do_not_dwnl_cache,
        load=False,
        install_dir=args.install_dir,
        logger=logger,
    )

    # try:
    if args.build:
        cache.Build(interactive=args.interactive)
    elif args.reaction_rules is not None:
        print_attr(
            cache, "rr_reactions", args.reaction_rules, args.do_not_dwnl_cache, logger
        )
    elif args.reactions is not None:
        print_attr(
            cache, "template_reactions", args.reactions, args.do_not_dwnl_cache, logger
        )
    elif args.compounds is not None:
        print_attr(cache, "cid_strc", args.compounds, args.do_not_dwnl_cache, logger)
    elif args.chem_xref is not None:
        print_attr(cache, "cid_xref", args.chem_xref, args.do_not_dwnl_cache, logger)
    else:
        cache.Load(
            interactive=args.interactive, do_not_dwnl_cache=args.do_not_dwnl_cache
        )


def print_attr(
    cache: "rrCache",
    attr: str,
    attr_lst: List,
    do_not_dwnl_cache: bool,
    logger: Logger = getLogger(__file__),
) -> None:
    cache.Load(attrs=[attr], do_not_dwnl_cache=do_not_dwnl_cache)
    if attr_lst == []:
        print(dumps(cache.get(attr), indent=4))
    else:
        for id in attr_lst:
            try:
                print(id + ":", dumps(cache.get(attr)[id], indent=4))
            except KeyError:
                logger.error(
                    "ID not found in rrCache('{attr}'): {id}".format(attr=attr, id=id)
                )


if __name__ == "__main__":
    entry_point()
