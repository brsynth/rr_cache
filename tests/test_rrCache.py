"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

# import logging
from unittest import TestCase
from os import (
    remove as os_rm,
    path as os_path
)
from logging import Logger
from json import load as json_load
from rr_cache import rrCache
from brs_utils import (
    create_logger,
    extract_gz,
    check_file_size
)

HERE = os_path.dirname(os_path.abspath(__file__))
DATA_PATH = os_path.join(HERE, 'data')


class Test_rrCache(TestCase):

    data_type = 'mnx3.1'
    outdir = f'cache-{data_type}'
    cache = rrCache(attrs=None, data_type=data_type, interactive=False)

    # Not possible to compare hashes since
    # files contain dict that have to be sorted
    # before comparing them and then fill up the memory
    # Size of gunzipped files

    def setUp(self, logger: Logger = None):
        if logger is None:
            self.logger = create_logger(__name__, 'ERROR')
        else:
            self.logger = logger
        # Set attributes from data files
        for elem in ['reactions', 'retrorules']:
            with open(os_path.join(DATA_PATH, f'{elem}.json'), 'r') as f:
                setattr(self, f'{elem}', json_load(f))
        for elem in ['compounds', 'metrics']:
            with open(os_path.join(DATA_PATH, f'{elem}_{self.data_type}.json'), 'r') as f:
                setattr(self, f'{elem}', json_load(f))

    def test_all_attr(self):
        r"""Test of loading all attributes in rrCache and store them in files.

        Method: Load a full rrCache in 'file' store mode. Then, for each
        attribute, compare its length with it is supposed to be.
        """
        cache = rrCache(data_type=self.data_type, interactive=False, logger=self.logger)
        for attr in self.metrics:
            length = self.metrics[attr]['length']
            with self.subTest(attr=attr, length=length):
                self.assertEqual(len(cache.get(attr)), length)

    def test_single_attr_file(self):
        r"""Test of loading each attribute in rrCache and store it in a file.

        Method: Load a rrCache in 'file' store mode for each single attribute.
        Then, compare its length with it is supposed to be.
        """
        for attr in self.metrics:
            length = self.metrics[attr]['length']
            with self.subTest(attr=attr, length=length):
                cache = rrCache([attr], data_type=self.data_type, interactive=False, logger=self.logger)
                self.assertEqual(len(cache.get(attr)), length)

    def test_generate_cache(self):
        r"""Test of genrating all rrCache files from input_cache.

        Method: Generate a full rrCache. Then, for each file, compare its size
        with it is supposed to be.
        """
        self.skipTest("Too long, not in standard tests")
        rrCache.generate_cache(self.outdir, interactive=False, logger=self.logger)
        for name in self.metrics:
            filepath = os_path.join(self.outdir, f'{name}.json.gz')
            outfile = extract_gz(filepath, self.outdir)
            self.assertTrue(
                check_file_size(
                    outfile,
                    self.metrics[name]['file_size'],
                    self.logger
                )
            )
            os_rm(outfile)

    def test_get_compound(self):
        self.assertDictEqual(
            self.cache.get_compound('MNXM2'),
            self.compounds['MNXM2']
        )

    def test_get_list_of_compounds(self):
        self.assertTrue(
            'MNXM2' in self.cache.get_list_of_compounds()
        )
        self.assertEqual(
            len(self.cache.get_list_of_compounds()),
            self.metrics['cid_strc']['length']
        )

    def test_get_reaction(self):
        self.assertDictEqual(
            self.cache.get_reaction('MNXR94688'),
            self.reactions['MNXR94688']
        )
        self.assertDictEqual(
            self.cache.get_reaction('MNXR94688'),
            self.reactions['MNXR94688']
        )

    def test_get_list_of_reactions(self):
        self.assertTrue(
            'MNXR94688' in self.cache.get_list_of_reactions()
        )
        self.assertEqual(
            len(self.cache.get_list_of_reactions()),
            self.metrics['template_reactions']['length']
        )

        self.assertTrue(
            'MNXR94688' in self.cache.get_list_of_reactions()
        )
        self.assertEqual(
            len(self.cache.get_list_of_reactions()),
            self.metrics['template_reactions']['length']
        )

    def test_get_reaction_rule(self):
        self.assertDictEqual(
            self.cache.get_reaction_rule('RR-02-f85f00f767901186-16-F'),
            self.retrorules['RR-02-f85f00f767901186-16-F']
        )

    def test_get_list_of_reaction_rules(self):
        self.assertTrue(
            'RR-02-f85f00f767901186-16-F' in self.cache.get_list_of_reaction_rules()
        )
        self.assertEqual(
            len(self.cache.get_list_of_reaction_rules()),
            self.metrics['rr_reactions']['length']
        )
