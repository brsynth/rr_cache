"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

# import logging
from unittest import TestCase
from os import remove as os_rm, path as os_path
from logging import Logger
from json import load as json_load
from rr_cache import rrCache
from brs_utils import create_logger, extract_gz, check_file_size

HERE = os_path.dirname(os_path.abspath(__file__))
DATA_PATH = os_path.join(HERE, "data")


class Test_rrCache(TestCase):

    cspaces = ["mnx3.1"]

    # Not possible to compare hashes since
    # files contain dict that have to be sorted
    # before comparing them and then fill up the memory
    # Size of gunzipped files

    def setUp(self, logger: Logger = None):
        if logger is None:
            self.logger = create_logger(__name__, "ERROR")
        else:
            self.logger = logger
        # Set attributes from data files
        for elem in ["reactions", "retrorules"]:
            with open(os_path.join(DATA_PATH, f"{elem}.json"), "r") as f:
                setattr(self, f"{elem}", json_load(f))
        for elem in ["compounds", "metrics"]:
            for cspace in self.cspaces:
                with open(os_path.join(DATA_PATH, f"{elem}_{cspace}.json"), "r") as f:
                    setattr(self, f"{elem}_{cspace}", json_load(f))
        # Set caches
        for cspace in self.cspaces:
            cache = rrCache(cspace=cspace, interactive=False, logger=self.logger)
            setattr(self, f"cache_{cspace}", cache)

    def test_all_attr(self):
        r"""Test of loading all attributes in rrCache and store them in files.

        Method: Load a full rrCache in 'file' store mode. Then, for each
        attribute, compare its length with it is supposed to be.
        """
        for cspace in self.cspaces:
            cache = getattr(self, f"cache_{cspace}")
            metrics = getattr(self, f"metrics_{cspace}")
            for attr in metrics:
                length = metrics[attr]["length"]
                with self.subTest(attr=attr, length=length):
                    self.assertEqual(len(cache.get(attr)), length)

    def test_generate_cache(self):
        r"""Test of genrating all rrCache files from input_cache.

        Method: Generate a full rrCache. Then, for each file, compare its size
        with it is supposed to be.
        """
        self.skipTest("Too long, not in standard tests")
        for cspace in self.cspaces:
            outdir = f"cache-{cspace}"
            rrCache.generate_cache(outdir, interactive=False, logger=self.logger)
            for name in self.metrics:
                filepath = os_path.join(outdir, f"{name}.json.gz")
                outfile = extract_gz(filepath, outdir)
                self.assertTrue(
                    check_file_size(
                        outfile, self.metrics[name]["file_size"], self.logger
                    )
                )
                os_rm(outfile)

    def test_get_compound(self):
        for cspace in self.cspaces:
            cache = getattr(self, f"cache_{cspace}")
            compounds = getattr(self, f"compounds_{cspace}")
            self.assertDictEqual(cache.get_compound("MNXM2"), compounds["MNXM2"])

    def test_get_list_of_compounds(self):
        for cspace in self.cspaces:
            cache = getattr(self, f"cache_{cspace}")
            metrics = getattr(self, f"metrics_{cspace}")
            self.assertTrue("MNXM2" in cache.get_list_of_compounds())
            self.assertEqual(
                len(cache.get_list_of_compounds()), metrics["cid_strc"]["length"]
            )

    def test_get_reaction(self):
        for cspace in self.cspaces:
            cache = getattr(self, f"cache_{cspace}")
            self.assertDictEqual(
                cache.get_reaction("MNXR94688"), self.reactions["MNXR94688"]
            )
            self.assertDictEqual(
                cache.get_reaction("MNXR94688"), self.reactions["MNXR94688"]
            )

    def test_get_list_of_reactions(self):
        for cspace in self.cspaces:
            cache = getattr(self, f"cache_{cspace}")
            metrics = getattr(self, f"metrics_{cspace}")
            self.assertTrue("MNXR94688" in cache.get_list_of_reactions())
            self.assertEqual(
                len(cache.get_list_of_reactions()),
                metrics["template_reactions"]["length"],
            )

            self.assertTrue("MNXR94688" in cache.get_list_of_reactions())
            self.assertEqual(
                len(cache.get_list_of_reactions()),
                metrics["template_reactions"]["length"],
            )

    def test_get_reaction_rule(self):
        for cspace in self.cspaces:
            cache = getattr(self, f"cache_{cspace}")
            self.assertDictEqual(
                cache.get_reaction_rule("RR-02-f85f00f767901186-16-F"),
                self.retrorules["RR-02-f85f00f767901186-16-F"],
            )

    def test_get_list_of_reaction_rules(self):
        for cspace in self.cspaces:
            cache = getattr(self, f"cache_{cspace}")
            metrics = getattr(self, f"metrics_{cspace}")
            self.assertTrue(
                "RR-02-f85f00f767901186-16-F" in cache.get_list_of_reaction_rules()
            )
            self.assertEqual(
                len(cache.get_list_of_reaction_rules()),
                metrics["rr_reactions"]["length"],
            )
