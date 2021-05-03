"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

# import logging
from unittest  import TestCase
from rr_cache import rrCache
from brs_utils import (
    create_logger,
    extract_gz
)
from os import (
    remove as os_rm,
    path as os_path
)


class Test_rrCache(TestCase):


    def setUp(self):
        self.logger = create_logger(__name__, 'ERROR')

    def test_all_attr_db(self):
        r"""Test of loading all attributes in rrCache and store them in a db.

        Method: Load a full rrCache in 'db' store mode (localhost). Then, for
        each attribute, compare its length with it is supposed to be.
        """
        self.skipTest("Tool long, not in standard tests")
        cache = rrCache('localhost', logger=self.logger)
        for attr, length in self.attributes:
            with self.subTest(attr=attr, length=length):
                self.assertEqual(len(cache.get(attr)), length)

    def test_all_attr_file(self):
        r"""Test of loading all attributes in rrCache and store them in files.

        Method: Load a full rrCache in 'file' store mode. Then, for each
        attribute, compare its length with it is supposed to be.
        """
        cache = rrCache('file', logger=self.logger)
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
                cache = rrCache('file', [attr], logger=self.logger)
                self.assertEqual(len(cache.get(attr)), length)

    def test_generate_cache(self):
        r"""Test of genrating all rrCache files from input_cache.

        Method: Generate a full rrCache. Then, for each file, compare its size
        with it is supposed to be.
        """
        self.skipTest("Tool long, not in standard tests")
        cache.generate_cache(self.outdir, logger=self.logger)
        for file, size in self.files:
            outfile = extract_gz(file, self.outdir)
            self.assertTrue(Main._check_file_size(outfile, size))
            os_rm(outfile)


    outdir = 'cache-3.2'

    # Not possible to compare hashes since files contain dict that have to be sorted before comparing them and then fill up the memory
    # Size of gunzipped files
    files = [
    (os_path.join(outdir, 'chebi_cid.json.gz'              ), 2786801),
    (os_path.join(outdir, 'cid_name.json.gz'               ), 55787548),
    (os_path.join(outdir, 'cid_strc.json.gz'               ), 296896910),
    (os_path.join(outdir, 'cid_xref.json.gz'               ), 88383985),
    (os_path.join(outdir, 'comp_xref.json.gz'              ), 51059),
    (os_path.join(outdir, 'deprecatedCID_cid.json.gz'      ), 423443),
    (os_path.join(outdir, 'deprecatedCompID_compid.json.gz'), 89832),
    (os_path.join(outdir, 'deprecatedRID_rid.json.gz'      ), 1437122),
    (os_path.join(outdir, 'inchikey_cid.json.gz'           ), 20071352),
    (os_path.join(outdir, 'rr_full_reactions.json.gz'      ), 7643885),
    (os_path.join(outdir, 'rr_reactions.json.gz'           ), 84656878)
    ]

    attributes = [
    ('chebi_cid',               123835),
    ('cid_name',                691482),
    ('cid_strc',                655684),
    ('cid_xref',                691494),
    ('comp_xref',               40),
    ('deprecatedCID_cid',       16267),
    ('deprecatedCompID_compid', 4370),
    ('deprecatedRID_rid',       53413),
    ('inchikey_cid',            332146),
    ('rr_full_reactions',       41970),
    ('rr_reactions',            229862)
    ]
