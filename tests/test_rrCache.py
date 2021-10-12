"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

# import logging
from unittest  import TestCase
from rr_cache import rrCache
from brs_utils import (
    create_logger,
    extract_gz,
    check_file_size
)
from os import (
    remove as os_rm,
    path as os_path
)
from logging import Logger


class Test_rrCache(TestCase):

    def setUp(self, logger: Logger = None):
        if logger is None:
            self.logger = create_logger(__name__, 'ERROR')
        else:
            self.logger = logger

    def test_all_attr(self):
        r"""Test of loading all attributes in rrCache and store them in files.

        Method: Load a full rrCache in 'file' store mode. Then, for each
        attribute, compare its length with it is supposed to be.
        """
        cache = rrCache(logger=self.logger)
        for attr, length in self.attributes:
            with self.subTest(attr=attr, length=length):
                self.assertEqual(len(cache.get(attr)), length)

    def test_single_attr_file(self):
        r"""Test of loading each attribute in rrCache and store it in a file.

        Method: Load a rrCache in 'file' store mode for each single attribute.
        Then, compare its length with it is supposed to be.
        """
        for attr, length in self.attributes:
            with self.subTest(attr=attr, length=length):
                cache = rrCache([attr], logger=self.logger)
                self.assertEqual(len(cache.get(attr)), length)

    def test_generate_cache(self):
        r"""Test of genrating all rrCache files from input_cache.

        Method: Generate a full rrCache. Then, for each file, compare its size
        with it is supposed to be.
        """
        self.skipTest("Too long, not in standard tests")
        rrCache.generate_cache(self.outdir, logger=self.logger)
        for file, size in self.files:
            outfile = extract_gz(file, self.outdir)
            self.assertTrue(check_file_size(outfile, size, self.logger))
            os_rm(outfile)

    def test_get_compound(self):
        cache = rrCache(attrs=None)
        self.assertDictEqual(
            cache.get_compound('MNXM2'),
            self.__MNXM2
        )

    def test_get_list_of_compounds(self):
        cache = rrCache(attrs=None)
        self.assertTrue(
            'MNXM2' in cache.get_list_of_compounds()
        )
        self.assertEqual(
            len(cache.get_list_of_compounds()),
            655684
        )

    def test_get_reaction(self):
        cache = rrCache(attrs=None)
        self.assertDictEqual(
            cache.get_reaction('MNXR94688'),
            self.__MNXR94688
        )

    def test_get_list_of_reactions(self):
        cache = rrCache(attrs=None)
        self.assertTrue(
            'MNXR94688' in cache.get_list_of_reactions()
        )
        self.assertEqual(
            len(cache.get_list_of_reactions()),
            44045
        )

    def test_get_reaction_rule(self):
        cache = rrCache(attrs=None)
        self.assertDictEqual(
            cache.get_reaction_rule('RR-02-f85f00f767901186-16-F'),
            self.__RR_02_f85f00f767901186_16_F
        )

    def test_get_list_of_reaction_rules(self):
        cache = rrCache(attrs=None)
        self.assertTrue(
            'RR-02-f85f00f767901186-16-F' in cache.get_list_of_reaction_rules()
        )
        self.assertEqual(
            len(cache.get_list_of_reaction_rules()),
            229862
        )

    outdir = 'cache-3.2'

    # Not possible to compare hashes since files contain dict that have to be sorted before comparing them and then fill up the memory
    # Size of gunzipped files
    files = [
    (os_path.join(outdir, 'chebi_cid.json.gz'), 2786801),
    (os_path.join(outdir, 'cid_name.json.gz'), 55787548),
    (os_path.join(outdir, 'cid_strc.json.gz'), 296896910),
    (os_path.join(outdir, 'cid_xref.json.gz'), 88383985),
    (os_path.join(outdir, 'comp_xref.json.gz'), 51059),
    (os_path.join(outdir, 'deprecatedCID_cid.json.gz'), 423443),
    (os_path.join(outdir, 'deprecatedCompID_compid.json.gz'), 89832),
    (os_path.join(outdir, 'deprecatedRID_rid.json.gz'), 1437122),
    (os_path.join(outdir, 'inchikey_cid.json.gz'), 20071352),
    (os_path.join(outdir, 'template_reactions.json.gz'), 7643885),
    (os_path.join(outdir, 'rr_reactions.json.gz'), 84656878)
    ]

    attributes = [
    ('chebi_cid', 123835),
    ('cid_name', 691481),
    ('cid_strc', 655684),
    ('cid_xref', 691494),
    ('comp_xref', 40),
    ('deprecatedCID_cid', 16267),
    ('deprecatedCompID_compid', 4370),
    ('deprecatedRID_rid', 53413),
    ('inchikey_cid', 332146),
    ('template_reactions', 44045),
    ('rr_reactions', 229862)
    ]

    __MNXM2 = {
        'formula': 'H2O',
        'smiles': 'O',
        'inchi': 'InChI=1S/H2O/h1H2',
        'inchikey': 'XLYOFNOQVPJJNP-UHFFFAOYSA-N',
        'cid': 'MNXM2',
        'name': 'H2O'
    }

    __MNXR94688 = {
        'left': {
            'MNXM1': 1,
            'MNXM17': 1,
            'MNXM91349': 1
        },
        'right': {
            'MNXM52': 1
        },
        'direction': 0,
        'main_left': [
            'MNXM91349'
        ],
        'main_right': [
            'MNXM52'
        ]
    }

    __RR_02_f85f00f767901186_16_F = {
        'MNXR100137': {
            'rule_id': 'RR-02-f85f00f767901186-16-F',
            'rule_score': 0.7486085128675456,
            'reac_id': 'MNXR100137',
            'subs_id': 'MNXM2210',
            'rel_direction': 1,
            'left': {
                'MNXM2210': 1
            },
            'right': {
                'MNXM139': 1,
                'MNXM83': 1
            }
        },
        'MNXR112247': {
            'rule_id': 'RR-02-f85f00f767901186-16-F',
            'rule_score': 0.7486085128675456,
            'reac_id': 'MNXR112247',
            'subs_id': 'MNXM7283',
            'rel_direction': 1,
            'left': {
                'MNXM7283': 1
            },
            'right': {
                'MNXM2473': 1,
                'MNXM83': 1
            }
        },
        'MNXR116653': {
            'rule_id': 'RR-02-f85f00f767901186-16-F',
            'rule_score': 0.7486085128675456,
            'reac_id': 'MNXR116653',
            'subs_id': 'MNXM13314',
            'rel_direction': 1,
            'left': {
                'MNXM13314': 1
            },
            'right': {
                'MNXM13541': 1,
                'MNXM83': 1
            }
        }
    }
