from os import (
    path as os_path,
    makedirs
)
from tempfile import NamedTemporaryFile
from rdkit.Chem import (
    MolFromSmiles,
    MolFromInchi,
    MolToSmiles,
    MolToInchi,
    MolToInchiKey,
)
from csv import (
    DictReader as csv_DictReader,
    reader as csv_reader
)
from json import (
    dump as json_dump,
    load as json_load
)
from gzip       import open as gzip_open
from re         import findall as re_findall
# from time       import time as time_time
from requests   import exceptions as r_exceptions
from hashlib    import sha512
from pathlib    import Path
from colored    import (
    attr as c_attr,
)
from logging import (
    Logger,
    getLogger,
    StreamHandler
)
from typing import (
    List,
    Tuple,
    Dict
)
from brs_utils  import (
    print_start,
    print_progress,
    print_end,
    download,
    check_sha
)

from .Args import DEFAULTS


HERE = os_path.dirname(os_path.abspath( __file__ ))
DATA_PATH = os_path.join(HERE, 'data')


class FileCorruptedError(Exception):
    pass


class rrCache:
    """Class to generate the cache

    Contains all the functions that parse different files, used to calculate the thermodynamics and the FBA of the the other steps. These should be called only when the files have changes

    """

    # __input__cache_url = 'ftp://ftp.vital-it.ch/databases/metanetx/MNXref/3.2/'

    # static attribues
    with open(os_path.join(DATA_PATH, 'convert.json'), 'r') as f:
        __convertMNXM = json_load(f)

    # name: sha512sum
    with open(os_path.join(DATA_PATH, 'input_cache.json'), 'r') as f:
        __input__cache = json_load(f)

    # Attributes with dependencies (other attributes + input_cache files)
    with open(os_path.join(DATA_PATH, 'cache.json'), 'r') as f:
        __cache = json_load(f)
    __attributes_list = list(__cache.keys())

    ## Cache constructor
    def __init__(
        self,
        attrs: List = DEFAULTS['attrs'],
        cache_dir: str = DEFAULTS['cache_dir'],
        mnx_version: str = DEFAULTS['mnx_version'],
        logger: Logger = getLogger(__name__)
    ) -> 'rrCache':

        self.logger = logger
        self.logger.debug('New instance of rrCache')
        self.logger.debug('attrs: '+str(attrs))
        self.logger.debug('cache_dir: '+str(cache_dir))
        self.logger.debug('mnx_version: '+str(mnx_version))

        self.__mnx_version = mnx_version
        if cache_dir == DEFAULTS['cache_dir']:
            self.__cache_dir = os_path.join(
                HERE,
                'cache',
                f'mnx_{self.__mnx_version}'
            )
        else:
            self.__cache_dir = cache_dir
        self.load(attrs)


    def load(self, attrs: List = DEFAULTS['attrs']):

        if attrs is None:
            return

        if attrs != []:
            if not isinstance(attrs, list):
                self.logger.warning('\'attrs\' argument is not of type list, trying to convert...')
                self.__attributes_list = [attrs]
                self.logger.warning('OK')
            else:
                self.__attributes_list = attrs

        for attr in self.__attributes_list:
            setattr(self, '__'+attr, None)

        rrCache._check_or_download_cache_to_disk(
            self.__cache_dir,
            self.__attributes_list,
            self.logger
        )

        try:
            self._check_or_load_cache()
        except (r_exceptions.RequestException,
                r_exceptions.InvalidSchema,
                r_exceptions.ConnectionError):
            rrCache.generate_cache(self.__cache_dir)
            self._check_or_load_cache()


    @staticmethod
    def _check_or_download_cache_to_disk(
        cache_dir: str,
        attributes_list: Dict,
        logger: Logger = getLogger(__name__)
    ) -> None:
        logger.debug('cache_dir: '+str(cache_dir))
        logger.debug('attributes: '+str(attributes_list))

        print_start(logger, 'Downloading cache')

        for attr in attributes_list:
            print_progress()
            filename = rrCache.__cache[attr]['file']['name']
            full_filename = os_path.join(cache_dir, filename)

            try:
                if os_path.exists(full_filename):
                    fingerprint = rrCache.__cache[attr]['file']['fingerprint']
                    if check_sha(
                        full_filename,
                        fingerprint
                    ):
                        logger.debug(filename+" already downloaded")
                    else:  # sha not ok
                        logger.debug(
                            '\nfilename: ' + filename
                        + '\nlocation: ' + cache_dir
                        + '\nsha (computed): ' + sha512(Path(full_filename).read_bytes()).hexdigest()
                        + '\nsha (expected): ' + fingerprint
                        )
                        raise FileNotFoundError
                else:
                    raise FileNotFoundError

            except FileNotFoundError:
                logger.debug("Downloading "+filename+"...")
                # start_time = time_time()
                if not os_path.isdir(cache_dir):
                    makedirs(cache_dir, exist_ok=True)
                download(
                    rrCache.__cache[attr]['file']['url']+filename,
                    full_filename
                )
                # rrCache.__cache[attr] = True
                # end_time = time_time()

        print_end(logger)

    def get(self, attr: str):
        self.logger.debug(f'Getting attribute: {attr}')
        try:
            return getattr(self, '__'+attr)
        except Exception as e:
            self.logger.error(str(e))
            return None

    def __hasattr(self, attr: str):
        return hasattr(self, '__'+attr)

    def get_compound(self, cid: str):
        return self.__get_object('cid_strc', cid)

    def get_list_of_compounds(self):
        return self.__get_list_of_objects('cid_strc')

    def get_reaction(self, rxn_id: str):
        return self.__get_object('template_reactions', rxn_id)

    def get_list_of_reactions(self):
        return self.__get_list_of_objects('template_reactions')

    def get_reaction_rule(self, rr_id: str):
        return self.__get_object('rr_reactions', rr_id)

    def get_list_of_reaction_rules(self):
        return self.__get_list_of_objects('rr_reactions')

    def __get_object(self, attr: str, id: str):
        try:
            if not self.__hasattr(attr):
                self.load(attrs=[attr])
            return self.get(attr)[id]
        except Exception as e:
            self.logger.error(str(e))

    def __get_list_of_objects(self, attr: str):
        try:
            if not self.__hasattr(attr):
                self.load(attrs=[attr])
            return self.get(attr).keys()
        except Exception as e:
            self.logger.error(str(e))

    def set(self, attr: str, val: object):
        try:
            return setattr(self, '__'+attr, val)
        except Exception as e:
            self.logger.error(str(e))

    #####################################################
    ################# ERROR functions ###################
    #####################################################

    class Error(Exception):
        """Error function for the convertion of structures
        """
        pass


    class DepictionError(Error):
        """Error function for the convertion of structures
        """
        def __init__(self, message):
            """Constructor for the class

            :param message: The error handling message string

            :type message: str

            :rtype: None
            :return: None
            """
            # self.expression = expression
            self.message = message

    @staticmethod
    def generate_cache(
        outdir: str = DEFAULTS['cache_dir'],
        mnx_version: str = DEFAULTS['mnx_version'],
        logger: Logger = getLogger(__name__)
    ) -> None:

        print(f'Generating cache using MetaNetX version {mnx_version}')

        if outdir is DEFAULTS['cache_dir']:
            outdir = HERE
        input_cache_dir = os_path.join(
            outdir,
            'input-cache',
            f'mnx_{mnx_version}'
        )
        cache_dir = os_path.join(
            outdir,
            'cache',
            f'mnx_{mnx_version}'
        )
        if not os_path.isdir(cache_dir):
            makedirs(cache_dir)
        if not os_path.isdir(input_cache_dir):
            makedirs(input_cache_dir)

        # FETCH INPUT_CACHE FILES
        print_start(logger, 'Downloading input cache')
        for input_type, input in rrCache.__input__cache.items():
            # ignore MNX versions other that specified
            if input_type.startswith('mnx_') and \
            input_type != f'mnx_{mnx_version}':
                pass
            else:
                logger.debug(f'Downloading {input_type}...')
                for filename, fingerprint in input['files'].items():
                    rrCache._download_input_cache(
                        url=input['url'],
                        file=filename,
                        outdir=input_cache_dir,
                        fingerprint=fingerprint,
                        logger=logger
                    )

        print_end(logger)

        # GENERATE CACHE FILES AND STORE THEM TO DISK
        print_start(logger, 'Generating cache')
        deprecatedCID_cid = rrCache._gen_deprecatedCID_cid(input_cache_dir, cache_dir, logger)
        print_progress(logger)
        cid_strc, cid_name = rrCache._gen_cid_strc_cid_name(input_cache_dir, cache_dir, deprecatedCID_cid, logger)
        print_progress(logger)
        rrCache._gen_inchikey_cid(input_cache_dir, cache_dir, cid_strc, logger)
        print_progress(logger)
        del cid_strc, cid_name
        cid_xref = rrCache._gen_cid_xref(input_cache_dir, cache_dir, deprecatedCID_cid, logger)
        print_progress(logger)
        rrCache._gen_chebi_cid(input_cache_dir, cache_dir, cid_xref)
        print_progress(logger)
        del cid_xref
        deprecatedRID_rid = rrCache._gen_deprecatedRID_rid(input_cache_dir, cache_dir, logger)
        print_progress(logger)
        rrCache._gen_rr_reactions(input_cache_dir, cache_dir, logger)  # , deprecatedCID_cid, deprecatedRID_rid, logger)
        print_progress(logger)
        rrCache._gen_comp_xref_deprecatedCompID_compid(input_cache_dir, cache_dir, logger)
        print_progress(logger)
        rrCache._gen_template_reactions(input_cache_dir, cache_dir, deprecatedRID_rid, logger)  # , deprecatedCID_cid, deprecatedRID_rid, logger)
        print_progress(logger)
        del deprecatedCID_cid, deprecatedRID_rid
        print_progress(logger)
        print_end(logger)


    @staticmethod
    def _gen_deprecatedCID_cid(
        input_dir: str,
        outdir: str,
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        attribute = 'deprecatedCID_cid'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        deprecatedCID_cid = None
        f_deprecatedCID_cid = os_path.join(outdir, rrCache.__cache[attribute]['file']['name'])

        # Do not checksum since it is a dictionary
        if os_path.exists(f_deprecatedCID_cid) and check_sha(
            f_deprecatedCID_cid,
            rrCache.__cache[attribute]
        ):
            deprecatedCID_cid = rrCache._load_cache_from_file(f_deprecatedCID_cid)
            logger.debug("   Cache file already exists")
        else:
            logger.debug("   Generating data...")
            deprecatedCID_cid = rrCache._m_deprecatedMNXM(
                os_path.join(input_dir, 'chem_xref.tsv')
            )
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(deprecatedCID_cid, f_deprecatedCID_cid)

        return {
            'attr': deprecatedCID_cid,
            'file': f_deprecatedCID_cid
        }


    @staticmethod
    def _gen_cid_strc_cid_name(
        input_dir: str,
        outdir: str,
        deprecatedCID_cid: Dict,
        logger: Logger = getLogger(__name__)
    ) -> Dict:

        attribute = 'cid_strc, cid_name'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        cid_strc = None
        cid_name = None
        f_cid_strc = os_path.join(outdir, rrCache.__cache['cid_strc']['file']['name'])
        f_cid_name = os_path.join(outdir, rrCache.__cache['cid_name']['file']['name'])
        f_warnings = os_path.join(outdir, 'warnings.json.gz')

        # Do not checksum since it is a dictionary
        if os_path.exists(f_cid_strc) and check_sha(
            f_cid_strc,
            rrCache.__cache['cid_strc']
        ) and os_path.exists(f_cid_name) and check_sha(
            f_cid_name,
            rrCache.__cache['cid_name']
        ):
            cid_strc = rrCache._load_cache_from_file(f_cid_strc)
            logger.debug("   Cache file already exists")
        else:
            if not deprecatedCID_cid['attr']:
                logger.debug("   Loading input data from file...")
                deprecatedCID_cid = rrCache._load_cache_from_file(deprecatedCID_cid['file'])
                # print_OK()
            logger.debug("   Generating data...")
            cid_strc, cid_name, warnings = rrCache._m_mnxm_strc(
                os_path.join(input_dir, 'compounds.tsv.gz'),
                os_path.join(input_dir, 'chem_prop.tsv'),
                deprecatedCID_cid['attr']
            )
            # Replace compound IDs that have no structure with one that has.
            # Done from a manually built file
            with open(os_path.join(input_dir, 'MNXM_replacement_20190524.csv')) as csv_file:
                reader = csv_reader(csv_file, delimiter=' ')
                for row in reader:
                    if not row[0].startswith('#') and len(row) > 1:
                        if row[1] != 'R_group':
                            cid_strc[row[0]] = cid_strc[row[1]]
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(cid_strc, f_cid_strc)
            rrCache._store_cache_to_file(cid_name, f_cid_name)
            rrCache._store_cache_to_file(warnings, f_warnings)

        return {
            'attr': cid_strc,
            'file': f_cid_strc
        }, {
            'attr': cid_name,
            'file': f_cid_name
        }

    @staticmethod
    def _gen_inchikey_cid(
        input_dir: str,
        outdir: str,
        cid_strc: Dict,
        logger: Logger = getLogger(__name__)
    ) -> None:
        attribute = 'inchikey_cid'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        inchikey_cid = None
        f_inchikey_cid = os_path.join(outdir, rrCache.__cache[attribute]['file']['name'])

        # Do not checksum since it is a dictionary
        if os_path.exists(f_inchikey_cid) and check_sha(
            f_inchikey_cid,
            rrCache.__cache[attribute]
        ):
            logger.debug("   Cache file already exists")
        else:
            if not cid_strc['attr']:
                logger.debug("   Loading input data from file...")
                cid_strc['attr'] = rrCache._load_cache_from_file(cid_strc['file'])
            logger.debug("   Generating data...")
            inchikey_cid = rrCache._m_inchikey_cid(cid_strc['attr'])
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(inchikey_cid, f_inchikey_cid)
            del inchikey_cid


    @staticmethod
    def _gen_cid_xref(
        input_dir: str,
        outdir: str,
        deprecatedCID_cid: Dict,
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        attribute = 'cid_xref'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        cid_xref = None
        f_cid_xref = os_path.join(outdir, rrCache.__cache[attribute]['file']['name'])

        # Do not checksum since it is a dictionary
        if os_path.exists(f_cid_xref) and check_sha(
            f_cid_xref,
            rrCache.__cache[attribute]
        ):
            cid_xref = rrCache._load_cache_from_file(f_cid_xref)
            logger.debug("   Cache file already exists")
        else:
            if not deprecatedCID_cid['attr']:
                logger.debug("   Loading input data from file...")
                deprecatedCID_cid['attr'] = rrCache._load_cache_from_file(deprecatedCID_cid['file'])
            logger.debug("   Generating data...")
            cid_xref = rrCache._m_mnxm_xref(
                os_path.join(input_dir, 'chem_xref.tsv'),
                deprecatedCID_cid['attr']
            )
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(cid_xref, f_cid_xref)

        return {
            'attr': cid_xref,
            'file': f_cid_xref
        }


    @staticmethod
    def _gen_chebi_cid(
        input_dir: str,
        outdir: str,
        cid_xref: Dict,
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        attribute = 'chebi_cid'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        chebi_cid = None
        f_chebi_cid = os_path.join(outdir, rrCache.__cache[attribute]['file']['name'])

        # Do not checksum since it is a dictionary
        if os_path.exists(f_chebi_cid) and check_sha(
            f_chebi_cid,
            rrCache.__cache[attribute]
        ):
            logger.debug("   Cache file already exists")
        else:
            logger.debug("   Generating data...")
            chebi_cid = rrCache._m_chebi_cid(cid_xref['attr'])
            # print_OK()
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(chebi_cid, f_chebi_cid)
            del chebi_cid
            # print_OK()


    @staticmethod
    def _gen_deprecatedRID_rid(
        input_dir: str,
        outdir: str,
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        attribute = 'deprecatedRID_rid'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        deprecatedRID_rid = None
        f_deprecatedRID_rid = os_path.join(outdir, rrCache.__cache[attribute]['file']['name'])

        # Do not checksum since it is a dictionary
        if os_path.exists(f_deprecatedRID_rid) and check_sha(
            f_deprecatedRID_rid,
            rrCache.__cache[attribute]
        ):
            deprecatedRID_rid = rrCache._load_cache_from_file(f_deprecatedRID_rid)
            logger.debug("   Cache file already exists")
        else:
            logger.debug("   Generating data...")
            deprecatedRID_rid = rrCache._m_deprecatedMNXR(
                os_path.join(input_dir, 'reac_xref.tsv')
            )
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(deprecatedRID_rid, f_deprecatedRID_rid)

        return {
            'attr': deprecatedRID_rid,
            'file': f_deprecatedRID_rid
        }


    @staticmethod
    def _gen_rr_reactions(
        input_dir: str,
        outdir: str,
        # deprecatedCID_cid: Dict,
        # deprecatedRID_rid: Dict,
        logger: Logger = getLogger(__name__)
    ) -> None:
        attribute = 'rr_reactions'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        rr_reactions = None
        f_rr_reactions = os_path.join(outdir, rrCache.__cache[attribute]['file']['name'])

        # Do not checksum since it is a dictionary
        if os_path.exists(f_rr_reactions) and check_sha(
            f_rr_reactions,
            rrCache.__cache[attribute]
        ):
            logger.debug("   Cache file already exists")
        else:
            # if not deprecatedCID_cid['attr']:
            #     logger.debug("   Loading input data from file...")
            #     deprecatedCID_cid['attr'] = rrCache._load_cache_from_file(deprecatedCID_cid['file'])
            # if not deprecatedRID_rid['attr']:
            #     logger.debug("   Loading input data from file...")
            #     deprecatedRID_rid['attr'] = rrCache._load_cache_from_file(deprecatedRID_rid['file'])
            #     print_OK()
            logger.debug("   Generating data...")
            rr_reactions = rrCache._m_rr_reactions(
                os_path.join(input_dir, 'retrorules_rr02_flat_all.tsv.gz'),
                logger=logger
                # deprecatedCID_cid,
                # deprecatedRID_rid
            )
            # del deprecatedRID_rid
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(rr_reactions, f_rr_reactions)
            del rr_reactions


    @staticmethod
    def _gen_comp_xref_deprecatedCompID_compid(
        input_dir: str,
        outdir: str,
        logger: Logger = getLogger(__name__)
    ) -> None:
        attribute = 'comp_xref, deprecatedCompID_compid'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        comp_xref = deprecatedCompID_compid = None
        f_comp_xref = os_path.join(outdir, rrCache.__cache['comp_xref']['file']['name'])
        f_deprecatedCompID_compid = os_path.join(outdir, rrCache.__cache['deprecatedCompID_compid']['file']['name'])

        # Do not checksum since it is a dictionary
        if os_path.exists(f_comp_xref) and check_sha(
            f_comp_xref,
            rrCache.__cache['comp_xref']
        ) and os_path.exists(f_deprecatedCompID_compid) and check_sha(
            f_deprecatedCompID_compid,
            rrCache.__cache['deprecatedCompID_compid']
        ):
            logger.debug("   Cache files already exist")
            # print_OK()
        else:
            logger.debug("   Generating data...")
            comp_xref, deprecatedCompID_compid = rrCache._m_mnxc_xref(
                os_path.join(input_dir, 'comp_xref.tsv')
            )
            # print_OK()
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(comp_xref, f_comp_xref)
            # print_OK()
            del comp_xref
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(deprecatedCompID_compid, f_deprecatedCompID_compid)
            # print_OK()
            del deprecatedCompID_compid


    @staticmethod
    def _gen_template_reactions(
        input_dir: str,
        outdir: str,
        deprecatedRID_rid: Dict,
        logger: Logger = getLogger(__name__)
    ) -> None:
        logger.debug('Generating template_reactions')
        logger.debug('input_dir: '+str(input_dir))
        logger.debug('outdir: '+str(outdir))
        attribute = 'template_reactions'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        template_reactions = None
        f_template_reactions = os_path.join(outdir, rrCache.__cache[attribute]['file']['name'])

        # if os_path.exists(f_template_reactions) and check_sha(
        #     f_template_reactions,
        #     rrCache.__cache[attribute]
        # ):
        #     logger.debug("   Cache file already exists")
        # else:
        logger.debug("   Generating data...")
        template_reactions = rrCache._m_template_reactions(
            os_path.join(input_dir, 'rxn_recipes.tsv.gz')
        )
        for depRID, newRID in deprecatedRID_rid['attr'].items():
            try:
                template_reactions[depRID] = template_reactions[newRID]
            except KeyError as key:
                logger.warning(f'Reaction ID {key} not found in rxn_recipes.tsv.gz')
        logger.debug("   Writing data to file...")
        rrCache._store_cache_to_file(template_reactions, f_template_reactions)
        del template_reactions


    # def _load_from_file(self, filename):
    #     self.logger.debug("Loading "+filename+"...")
    #     data = self._load_cache_from_file(
    #         os_path.join(self.__cache_dir, filename)
    #     )
    #     return data


    def _check_or_load_cache(self):
        self._check_or_load_cache_in_memory()


    def _check_or_load_cache_in_memory(self):
        print_start(self.logger, 'Loading cache in memory')
        for attribute in self.__attributes_list:
            filename = rrCache.__cache[attribute]['file']['name']
            if self.get(attribute) is None:
                self.set(
                    attribute,
                    self._load_cache_from_file(
                        os_path.join(self.__cache_dir, filename)
                    )
                )
                print_progress(self.logger)
            else:
                self.logger.debug(attribute+" already loaded in memory")
        print_end(self.logger)


    @staticmethod
    def _download_input_cache(
        url: str,
        file: str,
        outdir: str,
        fingerprint: str=None,
        logger: Logger = getLogger(__name__)
    ):
        if not os_path.isdir(outdir):
            makedirs(outdir, exist_ok=True)
        filename = os_path.join(outdir, file)
        if (
            not os_path.isfile(filename)
            or not check_sha(filename, fingerprint)
        ):
            # start_time = time_time()
            rrCache.__download_input_cache(url, file, outdir)
            print_progress(logger)
            # end_time = time_time()
        if not check_sha(
            filename,
            fingerprint
        ):  # sha not ok
            logger.debug(f'\n\
                filename: {filename}\n\
                sha (computed): {sha512(Path(filename).read_bytes()).hexdigest()}\n\
                sha (expected): {fingerprint}\n\
            '
            )
            logger.error(f'\nUnable to download input-cache file {file}.')
            logger.error('\nEither the URL is broken or the file content has changed.')
            logger.error('\nExiting...\n')
            exit()
            raise FileCorruptedError


    @staticmethod
    def __download_input_cache(
        url: str,
        file: str,
        outdir: str,
        logger: Logger = getLogger(__name__)
    ):

        if not os_path.isdir(outdir):
            makedirs(outdir ,exist_ok=True)

        logger.debug(f'Downloading {file} from {url}')
        download(url+file, os_path.join(outdir, file))

    ##########################################################
    ################## Private Functions #####################
    ##########################################################

    ## Method to load data from file
    #
    #  Load data from file
    #
    #  @param self Object pointer
    #  @param filename File to fetch data from
    #  @return file content


    @staticmethod
    def _load_cache_from_file(filename, logger: Logger = getLogger(__name__)):
        logger.debug(filename)
        if filename.endswith('.gz') or filename.endswith('.zip'):
            fp = gzip_open(filename, 'rt', encoding='ascii')
        else:
            fp = open(filename, 'r')
        return json_load(fp)

    ## Method to store data into file
    #
    # Store data into file as json (to store dictionnary structure)
    #
    #  @param self Object pointer
    #  @param data Data to write into file
    #  @param filename File to write data into
    @staticmethod
    def _store_cache_to_file(data, filename):
        if filename.endswith('.gz') or filename.endswith('.zip'):
            fp = gzip_open(filename, 'wt', encoding='ascii')
        else:
            fp = open(filename, 'w')
        json_dump(data, fp)

    ## Function to create a dictionnary of old to new chemical id's
    #
    #  Generate a one-to-one dictionnary of old id's to new ones. Private function
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    @staticmethod
    def _checkCIDdeprecated(cid, deprecatedCID_cid):
        try:
            return deprecatedCID_cid[cid]
        except (KeyError, TypeError):
            return cid

    ## Function to create a dictionnary of old to new reaction id's
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    @staticmethod
    def _checkRIDdeprecated(rid, deprecatedRID_rid):
        try:
            return deprecatedRID_rid[rid]
        except (KeyError, TypeError):
            return rid

    #################################################################
    ################## Public functions #############################
    #################################################################


    ########################### MNX parsers #############################

    ## Function to parse the chem_xref.tsv and reac_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's.
    # This can include more than one old id per new one and thus returns a dictionnary. Private function
    #
    # @param xref_path Input file path
    # @return Dictionnary of identifiers
    # TODO: save the self.deprecatedCID_cid to be used in case there rp_paths uses an old version of MNX


    @staticmethod
    def _deprecatedMNX(xref_path):
        deprecatedMNX_mnx = {}
        with open(xref_path, 'rt') as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0] == '#':
                    mnx = row[0].split(':')
                    if mnx[0] == 'deprecated':
                        deprecatedMNX_mnx[mnx[1]] = row[1]
        return deprecatedMNX_mnx

    ## Status function that parses the chem_xref.tsv file for chemical cross-references and the different
    #
    # @param chem_xref_path Input file path
    # @return Dictionnary of chemical id to other chemical ids ex: deprecatedCID_cid['MNXM1'] = {'mnx': ['MNXM01', ...], ...}
    @staticmethod
    def _m_deprecatedMNXM(chem_xref_path):
        deprecatedCID_cid = {}
        deprecatedCID_cid = rrCache._deprecatedMNX(chem_xref_path)
        deprecatedCID_cid.update(rrCache.__convertMNXM)
        deprecatedCID_cid['MNXM01'] = 'MNXM1'
        return deprecatedCID_cid

    ## Function to parse the reac_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's.
    # This can include more than one old id per new one and thus returns a dictionnary. Private function
    #
    #  @param self Object pointer
    #  @param reac_xref_path Input file path
    #  @return Dictionnary of identifiers
    @staticmethod
    def _m_deprecatedMNXR(reac_xref_path):
        return rrCache._deprecatedMNX(reac_xref_path)

    ## Function to parse the chemp_prop.tsv file from MetanetX and compounds.tsv from RetroRules. Uses the InchIkey as key to the dictionnary
    #
    #  Generate a dictionnary gaving the formula, smiles, inchi and inchikey for the components
    # TODO: Seperate this function to parse the chem_prop (mnx specific) and the compounds.tsv from RetroRules (generic, not mnx specific)
    # Structure of return: cid_strc['MNXM1'] = {'formula': 'H', 'smiles': '[H+]', 'inchi': 'InChI=1S/p+1', 'inchikey': 'GPRLSGONYQIRFK-UHFFFAOYSA-N'}
    #
    #  @param rr_compounds_path Path to the RetroRules file
    #  @param chem_prop_path Path to the chem_prop.tsv file
    #  @param deprecatedCID_cid Dictionnary of deprecated CID to cid
    #  @return cid_strc Dictionnary of formula, smiles, inchi and inchikey
    @staticmethod
    def _m_mnxm_strc(
        rr_compounds_path: str,
        chem_prop_path: str,
        deprecatedCID_cid: Dict,
        logger: Logger = getLogger(__name__)
    ) -> Tuple[Dict, Dict]:

        cid_strc = {}
        cid_name = {}

        for row in csv_DictReader(gzip_open(rr_compounds_path, 'rt'), delimiter='\t'):
            tmp = {
                'formula':  None,
                'smiles':   None,
                'inchi':    row['inchi'],
                'inchikey': None,
                'cid':      rrCache._checkCIDdeprecated(row['cid'], deprecatedCID_cid),
                'name':     None
            }
            try:
                resConv = rrCache._convert_depiction(idepic=tmp['inchi'], itype='inchi', otype={'smiles', 'inchikey'})
                for i in resConv:
                    tmp[i] = resConv[i]
            except rrCache.DepictionError as e:
                logger.warning('Could not convert some of the structures: '+str(tmp))
                logger.warning(e)
            cid_strc[tmp['cid']] = tmp

        # logger.addHandler(
        #     FileHandler(
        #         os_path.join('warnings.log'), 'w'
        #     )
        # )
        with open(chem_prop_path, 'rt') as f:
            # read CSV with both tab and space as delimiters
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if row[0].startswith('#'):
                    header = row
                else:
                    for i in range(len(header)):
                        # remove '#' from column field and
                        # convert to lower case
                        field = header[i].replace('#', '').lower()
                        tmp[field] = row[i]
                    mnxm = rrCache._checkCIDdeprecated(row[0], deprecatedCID_cid)
                    # tmp = {
                    #     'formula':  row[2],
                    #     'smiles': row[6],
                    #     'inchi': row[5],
                    #     'inchikey': row[8],
                    #     'cid': mnxm,
                    #     'name': row[1]
                    # }
                    for i in tmp:
                        if tmp[i] == '' or tmp[i] == 'NA':
                            tmp[i] = None
                    if mnxm not in cid_name and tmp['name']:
                        cid_name[mnxm] = tmp['name']
                    if mnxm in cid_strc:
                        cid_strc[mnxm]['formula'] = row[2]
                        cid_strc[mnxm]['name'] = row[1]
                        if not cid_strc[mnxm]['smiles'] and tmp['smiles']:
                            cid_strc[mnxm]['smiles'] = tmp['smiles']
                        if not cid_strc[mnxm]['inchikey'] and tmp['inchikey']:
                            cid_strc[mnxm]['inchikey'] = tmp['inchikey']
                    else:
                        # check to see if the inchikey is valid or not
                        otype = set({})
                        if not tmp['inchikey']:
                            otype.add('inchikey')
                        if not tmp['smiles']:
                            otype.add('smiles')
                        if not tmp['inchi']:
                            otype.add('inchi')
                        itype = ''
                        # Check if the current compound has
                        # InChI or SMILES description
                        if tmp['inchi']:
                            itype = 'inchi'
                        elif tmp['smiles']:
                            itype = 'smiles'
                        else:
                            ter = StreamHandler.terminator
                            StreamHandler.terminator = "\n"
                            logger.warning('No InChI or SMILES for '+str(tmp))
                            StreamHandler.terminator = ter
                            continue
                        ter = StreamHandler.terminator
                        StreamHandler.terminator = "\n"
                        try:
                            resConv = rrCache._convert_depiction(idepic=tmp[itype], itype=itype, otype=otype)
                            for i in resConv:
                                tmp[i] = resConv[i]
                            logger.warning('Sructure conversion OK: '+str(tmp))
                        except rrCache.DepictionError as e:
                            logger.warning('Structure conversion FAILED: '+str(tmp))
                            logger.warning(e)
                        # StreamHandler.terminator = ter
                        cid_strc[tmp['cid']] = tmp
        # logger.removeHandler(logger.handlers[-1])

        return cid_strc, cid_name


    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of all cross references for a given chemical id (MNX) to other database id's
    #
    #  @param chem_xref_path MetaNetX chem_xref.tsv file path
    #  @param deprecatedCID_cid Dictionnary of deprecated chemical ids to uniform cid
    #  @return Dictionnary of cross references of a given chemical id
    # TODO: save the self.deprecatedCID_cid to be used in case there rp_paths uses an old version of MNX


    @staticmethod
    def _m_mnxm_xref(
        chem_xref_path: str,
        deprecatedCID_cid: Dict,
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        cid_xref = {}
        with open(chem_xref_path, 'rt') as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0] == '#':
                    mnx = rrCache._checkCIDdeprecated(row[1], deprecatedCID_cid)
                    if len(row[0].split(':')) == 1:
                        dbName = 'mnx'
                        dbId = row[0]
                    else:
                        dbName = row[0].split(':')[0]
                        dbId = ''.join(row[0].split(':')[1:])
                        if dbName == 'deprecated':
                            dbName = 'mnx'
                    # mnx
                    if mnx not in cid_xref:
                        cid_xref[mnx] = {}
                    if dbName not in cid_xref[mnx]:
                        cid_xref[mnx][dbName] = []
                    if dbId not in cid_xref[mnx][dbName]:
                        cid_xref[mnx][dbName].append(dbId)
                    ### DB ###
                    if dbName not in cid_xref:
                        cid_xref[dbName] = {}
                    if dbId not in cid_xref[dbName]:
                        cid_xref[dbName][dbId] = mnx
        return cid_xref


    ## Function to parse the comp_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of compartments id's (MNX) to other database id's
    #
    #  @param comp_xref_path The MetaNetX file that contains the cross references
    #  @return a The dictionnary of compartment identifiers
    # TODO: save the self.deprecatedCID_cid to be used in case there rp_paths uses an old version of MNX


    @staticmethod
    def _m_mnxc_xref(
        comp_xref_path,
        logger: Logger = getLogger(__name__)
    ) -> Tuple[Dict, Dict]:
        comp_xref = {}
        deprecatedCompID_compid = {}

        if not os_path.exists(comp_xref_path):
            logger.error('Could not read the file {filename}'.format(filename=comp_xref_path))
            return None

        with open(comp_xref_path, 'rt') as f:
            c = csv_reader(f, delimiter='\t')
            # not_recognised = []
            for row in c:
                # cid = row[0].split(':')
                if not row[0][0] == '#':
                    # collect the info
                    mnxc = row[1]
                    if len(row[0].split(':')) == 1:
                        dbName = 'mnx'
                        dbCompId = row[0]
                    else:
                        dbName = row[0].split(':')[0]
                        dbCompId = ''.join(row[0].split(':')[1:])
                        dbCompId = dbCompId.lower()
                    if dbName == 'deprecated':
                        dbName = 'mnx'
                    # create the dicts
                    if mnxc not in comp_xref:
                        comp_xref[mnxc] = {}
                    if dbName not in comp_xref[mnxc]:
                        comp_xref[mnxc][dbName] = []
                    if dbCompId not in comp_xref[mnxc][dbName]:
                        comp_xref[mnxc][dbName].append(dbCompId)
                    # create the reverse dict
                    if dbCompId not in deprecatedCompID_compid:
                        deprecatedCompID_compid[dbCompId] = mnxc

        return comp_xref, deprecatedCompID_compid


    ######################## RetroRules specific functions ##################


    ## Function to parse the rules_rall.tsv from RetroRules
    #
    #  Extract from the reactions rules the ruleID, the reactionID, the direction of the rule directed to the origin reaction
    #  Structure of the return: rr_reactions['RR-02-d2e7c5761b5a9b4b-04-F'] = {'MNXR139133': {'rule_id': 'RR-02-d2e7c5761b5a9b4b-04-F', 'rule_score': 0.3151075983206353, 'reac_id': 'MNXR139133', 'subs_id': 'MNXM89557', 'rel_direction': 1, 'left': {'MNXM89557': 1}, 'right': {'MNXM20': 1, 'MNXM722724': 1}}}
    #
    #  @param rules_rall_path Path to the RetroRules reaction rules
    #  @param deprecatedCID_cid Dictionnary of deprecated to uniformed chemical id's
    #  @param deprecatedRID_rid Dictionnary of deprecated to uniformed reaction id's
    #  @return Dictionnary describing each reaction rule


    @staticmethod
    def _m_rr_reactions(
        rules_rall_path: str,
        # deprecatedCID_cid,
        # deprecatedRID_rid,
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        rr_reactions = {}

        if not os_path.exists(rules_rall_path):
            logger.error('Could not read the rules_rall file ('+str(rules_rall_path)+')')
            return None

        for row in csv_DictReader(gzip_open(rules_rall_path, 'rt'), delimiter='\t'):
            # NOTE: as of now all the rules are generated using MNX
            # but it may be that other db are used, we are handling this case
            # WARNING: can have multiple products so need to seperate them
            products = {}

            for cid in row['Product_IDs'].split('.'):

                # cid = rrCache._checkCIDdeprecated(i, deprecatedCID_cid)
                if cid not in products:
                    products[cid] = 1
                else:
                    products[cid] += 1

            try:
                # WARNING: one reaction rule can have multiple reactions associated with them
                # To change when you can set subpaths from the mutliple numbers of
                # we assume that the reaction rule has multiple unique reactions associated
                if row['# Rule_ID'] not in rr_reactions:
                    rr_reactions[row['# Rule_ID']] = {}
                if row['# Rule_ID'] in rr_reactions[row['# Rule_ID']]:
                    logger.warning('There is already reaction '+str(row['# Rule_ID'])+' in reaction rule '+str(row['# Rule_ID']))
                rr_reactions[row['# Rule_ID']][row['Reaction_ID']] = {
                    'rule_id': row['# Rule_ID'],
                    'rule_score': float(row['Score_normalized']),
                    # 'reac_id': rrCache._checkRIDdeprecated(row['Reaction_ID'], deprecatedRID_rid),
                    # 'subs_id': rrCache._checkCIDdeprecated(row['Substrate_ID'], deprecatedCID_cid),
                    'reac_id': row['Reaction_ID'],
                    'subs_id': row['Substrate_ID'],
                    'rel_direction': int(row['Rule_relative_direction']),
                    # 'left': {rrCache._checkCIDdeprecated(row['Substrate_ID'], deprecatedCID_cid): 1},
                    'left': {row['Substrate_ID']: 1},
                    'right': products
                }

            except ValueError:
                logger.error('Problem converting rel_direction: '+str(row['Rule_relative_direction']))
                logger.error('Problem converting rule_score: '+str(row['Score_normalized']))

        return rr_reactions


    ## Generate complete reactions from the rxn_recipes.tsv from RetroRules
    #
    #  These are the compplete reactions from which the reaction rules are generated from. This is used to
    #  reconstruct the full reactions from monocomponent reactions
    #  Structure of the return: template_reactions['MNXR142257'] = {'left': {'MNXM4660': 1}, 'right': {'MNXM97172': 1}, 'direction': 0, 'main_left': ['MNXM4660'], 'main_right': ['MNXM97172']}
    #
    #  @param self The pointer object
    #  @param rxn_recipes_path Path to the recipes file
    #  @return Boolean that determines the success or failure of the function


    @staticmethod
    def _m_template_reactions(
        rxn_recipes_path: str,
        logger: Logger = getLogger(__name__)
    ) -> Dict:

        if not os_path.exists(rxn_recipes_path):
            logger.error('Cannot find file: '+str(rxn_recipes_path))
            return None

        reactions = {}

        for row in csv_DictReader(gzip_open(rxn_recipes_path, 'rt'), delimiter='\t'):

            # Read equation
            rxn = rrCache._read_equation(
                row['Equation'],
                row['#Reaction_ID'],
                logger
            )
            if rxn is None:
                # Pass to the next equation
                continue

            # Direction
            dir = rrCache._read_direction(
                row['Direction'],
                logger
            )
            if dir is None:
                # Pass to the next equation
                continue
            else:
                rxn['direction'] = dir

            # Others
            rxn['main_left'] = row['Main_left'].split(',')
            rxn['main_right'] = row['Main_right'].split(',')

            reactions[row['#Reaction_ID']] = rxn

        return reactions


    def _read_direction(
        dir: str,
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        try:
            _dir = int(dir)
        except ValueError:
            ter = StreamHandler.terminator
            StreamHandler.terminator = "\n"
            logger.warning(
                'Cannot convert direction value {dir} to int'.format(
                    dir=dir
                )
            )
            StreamHandler.terminator = ter
            # Pass to the next equation
            return None
        return _dir


    def _read_equation(
        eq: str,
        rxn_id: str,
        logger: Logger = getLogger(__name__)
    ) -> Dict:

        if not len(eq.split('=')) == 2:
            logger.warning('There should never be more or less than a left and right of an equation')
            logger.warning('Ignoring {eq}'.format(eq=eq))
            return None

        #### for character matching that are returned
        DEFAULT_STOICHIO_RESCUE = {
            '4n': 4, '3n': 3, '2n': 2, 'n': 1,
            '(n)': 1, '(N)': 1, '(2n)': 2, '(x)': 1,
            'N': 1, 'm': 1, 'q': 1,
            '0.01': 1, '0.1': 1, '0.5': 1, '1.5': 1,
            '0.02': 1, '0.2': 1,
            '(n-1)': 0, '(n-2)': -1
        }

        rxn = {}

        # Use int indices to be directly usable by re_findall
        # 0 = left, 1 = right
        for side in [0, 1]:
            rxn[side] = {}
            for spe in re_findall(r'(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+', eq.split('=')[side]):
                # 1) try to rescue if its one of the values
                try:
                    # rxn[side][rrCache._checkCIDdeprecated(spe[1], deprecatedCID_cid)] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                    rxn[side][spe[1]] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                except KeyError:
                    # 2) try to convert to int if its not
                    try:
                        # rxn[side][rrCache._checkCIDdeprecated(spe[1], deprecatedCID_cid)] = int(spe[0])
                        rxn[side][spe[1]] = float(spe[0])
                    except ValueError:
                        ter = StreamHandler.terminator
                        StreamHandler.terminator = "\n"
                        logger.warning(
                            f'Cannot convert stoichio coeff {spe[0]} in {rxn_id}'
                        )
                        StreamHandler.terminator = ter
                        # Stop parsing this equation and pass the next
                        return None

        return {
            'left': rxn[0],
            'right': rxn[1]
        }


    # ------------ GENERIC FUNCTIONS ------------ #

    # Convert chemical depiction to others type of depictions
    #
    # Usage example:
    # - convert_depiction(idepic='CCO', otype={'inchi', 'smiles', 'inchikey'})
    # - convert_depiction(
    #       idepic='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
    #       itype='inchi',
    #       otype={'inchi', 'smiles', 'inchikey'}
    #   )
    #
    #  @param self The object pointer
    #  @param idepic String depiction to be converted, str
    #  @param itype type of depiction provided as input, str
    #  @param otype types of depiction to be generated, {"", "", ..}
    #  @return odepic generated depictions, {"otype1": "odepic1", ..}
    @staticmethod
    def _convert_depiction(idepic, itype='smiles', otype={'inchikey'}):
        # Import (if needed)
        if itype == 'smiles':
            rdmol = MolFromSmiles(idepic, sanitize=True)
        elif itype == 'inchi':
            rdmol = MolFromInchi(idepic, sanitize=True)
        else:
            raise NotImplementedError(
                f'"{itype}" is not a valid input type'
            )
        if rdmol is None:  # Check imprt
            raise rrCache.DepictionError(
                f'Import error from depiction "{idepic}" of type "{itype}"'
            )
        # Export
        odepic = dict()
        for item in otype:
            if item == 'smiles':
                # MolToSmiles is tricky, one mays want to check the possible options..
                odepic[item] = MolToSmiles(rdmol)
            elif item == 'inchi':
                odepic[item] = MolToInchi(rdmol)
            elif item == 'inchikey':
                odepic[item] = MolToInchiKey(rdmol)
            else:
                raise NotImplementedError(f'"{otype}" is not a valid output type')
        return odepic


    # Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of all cross references
    #  for a given chemical id (MNX) to other database id's
    #  Structure if the return: chebi_cid['88281']: 'MXM2323'
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers
    # TODO: save the self.deprecatedCID_cid to be used in case
    #       there rp_paths uses an old version of MNX
    @staticmethod
    def _m_chebi_cid(cid_xref):
        chebi_cid = {}
        for cid in cid_xref:
            if 'chebi' in cid_xref[cid]:
                for c in cid_xref[cid]['chebi']:
                    chebi_cid[c] = cid
        return chebi_cid

    # Function to build the dictionnary to find the chemical id from inchikey
    #
    # @param cid_strc Dictionnary of chemical ID's to
    #        all the structure information associated with it
    # @return Dictionnary of InChIKey to chemical ID
    @staticmethod
    def _m_inchikey_cid(cid_strc):
        inchikey_cid = {}
        for cid in cid_strc:
            inchikey = cid_strc[cid]['inchikey']
            # This line is needed to put a value in 'inchikey',
            # otherwise there are some problems in future strucutres
            if not inchikey:
                inchikey = 'NO_INCHIKEY'
            if inchikey not in inchikey_cid:
                inchikey_cid[inchikey] = []
            inchikey_cid[inchikey].append(cid)
        return inchikey_cid
