from os import (
    path as os_path,
    rename as os_rename,
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
    dumps as json_dumps,
    load as json_load
)
from gzip import (
    open as gzip_open,
    GzipFile
)
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
from collections import (
    defaultdict,
    Counter
)
from brs_utils  import (
    print_start,
    print_progress,
    print_end,
    download,
    check_sha
)
from .Args import DEFAULTS, HERE, CONFIG_PATH


def ask_user_input():
    # Ask the user for input
    user_input = input("Press Enter to continue, 's' to stop, or 'c' to continue without prompting: ")
    # Check the user's input
    if user_input.lower() == 's':
        print("Stopping the program.")
        exit()
    elif user_input.lower() == 'c':
        print("Continuing without prompting.")
        return False
    else:
        print("Continuing...")
    return True


class FileCorruptedError(Exception):
    pass

class FingerprintError(Exception):
    pass


class rrCache:
    """Class to generate the cache

    Contains all the functions that parse different files, used to calculate the thermodynamics and the FBA of the the other steps. These should be called only when the files have changes

    """

    # __input__cache_url = 'ftp://ftp.vital-it.ch/databases/metanetx/MNXref/3.2/'


    ## Cache constructor
    def __init__(
        self,
        cspace: str = DEFAULTS['cspace'],
        interactive: bool = DEFAULTS['interactive'],
        do_not_dwnl_cache: bool = DEFAULTS['do_not_dwnl_cache'],
        load: bool = True,
        install_dir: str = DEFAULTS['install_dir'],
        logger: Logger = getLogger(__name__)
    ) -> 'rrCache':
        """Constructor for the class
        Args:
            cspace (str): Chemical space to use (e.g. mnx3.1, mnx4.4...).
            interactive (bool): Whether to ask the user for confirmation before overwriting existing files.
            do_not_dwnl_cache (bool): Whether to download the cache files from the internet.
            load (bool): Whether to load the cache files into memory.
            install_dir (str): Directory to install the cache.
            logger (Logger): Logger instance for logging messages.
        """

        self.logger = logger
        self.logger.debug('New instance of rrCache')
        self.logger.debug('cspace: '+str(cspace))
        self.logger.debug('interactive: '+str(interactive))
        self.logger.debug('do_not_dwnl_cache: '+str(do_not_dwnl_cache))
        self.logger.debug('load: '+str(load))
        self.logger.debug('install_dir: '+str(install_dir))

        self.__cspace = cspace

        # Config file
        # Search first in install_dir,
        # then in default config path
        cache_cfg_fln = os_path.join(install_dir, 'config', f'config_{self.__cspace}.json')
        if not os_path.exists(cache_cfg_fln):
            cache_cfg_fln = os_path.join(CONFIG_PATH, f'config_{self.__cspace}.json')

        try:
            with open(cache_cfg_fln, 'r') as f:
                cache_cfg = json_load(f)
                rrCache.__cache_sources = cache_cfg['sources']
                rrCache.__cache = cache_cfg['cache']
                rrCache.__type = cache_cfg.get('type', 'legacy')
        except FileNotFoundError:
            logger.error(f'Cache config file {cache_cfg_fln} not found, please check the --chemical-space argument')
            logger.error('Exiting...')
            exit(1)


        # Cache elements list
        rrCache.__attributes_list = list(rrCache.__cache.keys())

        # static attribues
        if 'static' in cache_cfg:
            try:
                convert_fln = os_path.join(CONFIG_PATH, cache_cfg['static'])
                with open(convert_fln, 'r') as f:
                    rrCache.__convertMNXM = json_load(f)
            except FileNotFoundError:
                rrCache.__convertMNXM = {}

        self.logger.info(f'Using {self.__cspace}')
        self.__input__cache_dir = os_path.join(install_dir, 'input-cache', self.__cspace)
        self.__cache_dir = os_path.join(install_dir, 'cache', self.__cspace)
        if load:
            self.Load(
                attrs=rrCache.__attributes_list,
                interactive=interactive,
                do_not_dwnl_cache=do_not_dwnl_cache,
                type=rrCache.__type
            )


    def get_cspace(self) -> str:
        """Get the chemical space used
        Returns:
            str: Chemical space
        """
        return self.__cspace
    
    def get_type(self) -> str:
        """Get the chemical space type used
        Returns:
            str: Chemical space type
        """
        return rrCache.__type

    def Load(
        self,
        attrs: List = [],
        interactive: bool = DEFAULTS['interactive'],
        do_not_dwnl_cache: bool = DEFAULTS['do_not_dwnl_cache'],
        type: str = 'legacy'
    ) -> None:
        """Load the cache attributes into memory
        Args:
            attrs (List): List of attributes to load, if None, all attributes are loaded
            interactive (bool): Whether to ask the user for confirmation before loading the cache
        Returns:
            None
        """
        self.logger.debug('Loading attributes: '+str(attrs))
        self.logger.debug('interactive: '+str(interactive))
        self.logger.debug('do_not_dwnl_cache: '+str(do_not_dwnl_cache))

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

        if not do_not_dwnl_cache:
            try:
                rrCache._check_or_download_cache_to_disk(
                    self.__cache_dir,
                    self.__attributes_list,
                    self.logger
                )
            except r_exceptions.MissingSchema as e:
                self.logger.warning(e)

        try:
            self._check_or_load_cache()
        except (FileNotFoundError,
                r_exceptions.RequestException,
                r_exceptions.InvalidSchema,
                r_exceptions.ConnectionError):
            self.Build(interactive=interactive)
            self._check_or_load_cache()


    def get_input_cache_dir(self) -> str:
        """Get the input cache directory
        Returns:
            str: Input cache directory
        """
        return self.__input__cache_dir

    @staticmethod
    def _check_or_download_cache_to_disk(
        cache_dir: str,
        attributes_list: Dict,
        logger: Logger = getLogger(__name__)
    ) -> None:
        logger.debug('cache_dir: '+str(cache_dir))
        logger.debug('attributes: '+str(attributes_list))
        # print_start(logger, 'Downloading cache')

        for attr in attributes_list:
            print_progress()
            filename = rrCache.__cache[attr]['file']['name']
            full_filename = os_path.join(cache_dir, filename)
            # try:
            if os_path.exists(full_filename):
                fingerprint = rrCache.__cache[attr]['file']['fingerprint']
                if check_sha(
                    full_filename,
                    fingerprint
                ):
                    logger.debug(filename+" already downloaded")
                else:  # sha not ok
                    logger.warning(
                        '\nfilename: ' + filename
                    + '\nlocation: ' + cache_dir
                    + '\nsha (computed): ' + sha512(Path(full_filename).read_bytes()).hexdigest()
                    + '\nsha (expected): ' + fingerprint
                    )
                    # raise FileNotFoundError
            else:
        #         raise FileNotFoundError

        # except FileNotFoundError:
                logger.debug("Downloading "+filename+"...")
                # start_time = time_time()
                if not os_path.isdir(cache_dir):
                    # cache_dir = '/mnx4.4'
                    makedirs(cache_dir, exist_ok=True)
                if rrCache.__cache[attr]['file']['url'] == "":
                    raise r_exceptions.MissingSchema(
                        f'URL for {attr} is empty, cannot download cache file'
                    )
                else:
                    download(
                        rrCache.__cache[attr]['file']['url']+filename,
                        full_filename
                    )
                # rrCache.__cache[attr] = True
                # end_time = time_time()

            # except FileNotFoundError:
            #     rrCache.__download_input_cache(rrCache.__cache[attr]['file']['url'], filename, cache_dir)
            #     print_progress(logger)

        # print_end(logger)

    def get(self, attr: str):
        # self.logger.debug(f'Getting attribute: {attr}')
        try:
            return getattr(self, '__'+attr)
        except AttributeError: # Try class attribute
            try:
                mangled = f"_{self.__class__.__name__}__{attr}"
                return getattr(self, mangled)
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
                self.Load(attrs=[attr])
            return self.get(attr)[id]
        except Exception as e:
            self.logger.error(str(e))

    def __get_list_of_objects(self, attr: str):
        try:
            if not self.__hasattr(attr):
                self.Load(attrs=[attr])
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

#    @staticmethod
    def Build(
        self,
        interactive: bool = DEFAULTS['interactive']
    ) -> None:
        """Generate the cache files and store them to disk.
        Args:
            interactive (bool): Whether to ask the user for confirmation before overwriting existing files.
            logger (Logger): Logger instance for logging messages.
        """
        self.logger.debug('interactive: '+str(interactive))

        # CACHE
        if os_path.isdir(self.__cache_dir):
            self.logger.warning(f'Cache directory {self.__cache_dir} already exists, data might be overwritten')
        else:
            makedirs(self.__cache_dir)

        # INPUT_CACHE
        if os_path.isdir(self.__input__cache_dir):
            self.logger.warning(f'Input cache directory {self.__input__cache_dir} already exists, existing data can be used, remove the directory to download new data')
        else:
            makedirs(self.__input__cache_dir)

        # FETCH INPUT_CACHE FILES
        print_start(self.logger, 'Checking input cache')
        for input_type, input in rrCache.__cache_sources.items():
            self.logger.debug(f'Checking {input_type}...')
            for filename, fingerprint in input['files'].items():
                # filename = 'metanetx' + os_path.sep + filename
                # Download if not exists or corrupted
                if not os_path.exists(os_path.join(self.__input__cache_dir, filename)):
                    self.logger.debug(f'{filename} not found in {self.__input__cache_dir}, downloading...')
                    rrCache._download_input_cache(
                        url=input['url'],
                        file=filename,
                        outdir=self.__input__cache_dir,
                        fingerprint=fingerprint,
                        logger=self.logger
                    )
                else:
                    if not check_sha(
                        os_path.join(self.__input__cache_dir, filename),
                        {
                            'file': {
                                'fingerprint': fingerprint
                            }
                        }
                    ):
                        self.logger.debug(f'{filename} found in input cache but corrupted, re-downloading...')
                        rrCache._download_input_cache(
                            url=input['url'],
                            file=filename,
                            outdir=self.__input__cache_dir,
                            fingerprint=fingerprint,
                            logger=self.logger
                        )
                    else:
                        self.logger.debug(f'{filename} found in input cache and valid, skipping download')
        print_end(self.logger)
        # BUILD CACHE FILES AND STORE THEM TO DISK
        print_start(self.logger, 'Building the cache')
        # try:
        #     deprecatedCID_cid = rrCache._gen_deprecatedCID_cid(self.__input__cache_dir, self.__cache_dir, self.logger)
        # except KeyError as e:
        #     self.logger.debug(f'{e} not found in input cache, skipping generation')
        #     deprecatedCID_cid = None
        # print_progress(self.logger)
        cid_strc, cid_name = rrCache._gen_cid_strc_cid_name(self.__input__cache_dir, self.__cache_dir, interactive=interactive, type=rrCache.__type, logger=self.logger)
        print_progress(self.logger)
        try:
            rrCache._gen_inchikey_cid(self.__input__cache_dir, self.__cache_dir, cid_strc, self.logger)
        except KeyError as e:
            self.logger.debug(f'{e} not found in input cache, skipping generation')
        print_progress(self.logger)
        del cid_strc, cid_name
        try:
            cid_xref = rrCache._gen_cid_xref(self.__input__cache_dir, self.__cache_dir, self.logger)
            print_progress(self.logger)
            rrCache._gen_chebi_cid(self.__input__cache_dir, self.__cache_dir, cid_xref)
            print_progress(self.logger)
            del cid_xref
        except KeyError as e:
            self.logger.debug(f'{e} not found in input cache, skipping generation')
        # try:
        #     deprecatedRID_rid = rrCache._gen_deprecatedRID_rid(self.__input__cache_dir, self.__cache_dir, self.logger)
        # except KeyError as e:
        #     deprecatedRID_rid = None
        #     self.logger.debug(f'{e} not found in input cache, skipping generation')
        # print_progress(self.logger)
        rrCache._gen_rr_reactions(self.__input__cache_dir, self.__cache_dir, type=rrCache.__type, logger=self.logger)  # , deprecatedCID_cid, deprecatedRID_rid, logger)
        print_progress(self.logger)
        try:
            rrCache._gen_comp_xref_deprecatedCompID_compid(self.__input__cache_dir, self.__cache_dir, self.logger)
        except KeyError as e:
            self.logger.debug(f'{e} not found in input cache, skipping generation')
        print_progress(self.logger)
        rrCache._gen_template_reactions(self.__input__cache_dir, self.__cache_dir, type=rrCache.__type, logger=self.logger)  # , deprecatedCID_cid, deprecatedRID_rid, logger)
        print_progress(self.logger)
        # del deprecatedCID_cid, deprecatedRID_rid
        print_progress(self.logger)
        print_end(self.logger)

    @staticmethod
    def _gen_deprecatedCID_cid(
        input_dir: str,
        outdir: str,
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        attribute = 'deprecatedCID_cid'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        logger.debug(f"   input_dir: {input_dir}")
        logger.debug(f"   outdir: {outdir}")
        deprecatedCID_cid = None
        f_deprecatedCID_cid = os_path.join(outdir, rrCache.__cache[attribute]['file']['name'])
        logger.debug(f"   f_deprecatedCID_cid: {f_deprecatedCID_cid}")

        # Do not checksum since it is a dictionary
        if os_path.exists(f_deprecatedCID_cid) and check_sha(
            f_deprecatedCID_cid,
            rrCache.__cache[attribute]
        ):
            deprecatedCID_cid = rrCache._load_json(f_deprecatedCID_cid)
            logger.debug("   Cache file already exists")
        else:
            logger.debug("   Generating data...")
            deprecatedCID_cid = rrCache._m_deprecatedMNXM(
                os_path.join(input_dir, 'chem_xref.tsv')
            )
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(deprecatedCID_cid, f_deprecatedCID_cid, logger=logger)

        return {
            'attr': deprecatedCID_cid,
            'file': f_deprecatedCID_cid
        }


    @staticmethod
    def _gen_cid_strc_cid_name(
        input_dir: str,
        outdir: str,
        interactive: bool = DEFAULTS['interactive'],
        type: str = 'legacy',
        logger: Logger = getLogger(__name__)
    ) -> Dict:

        attribute = 'cid_strc, cid_name'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        cid_strc = None
        # cid_name = None
        f_cid_strc = os_path.join(outdir, rrCache.__cache['cid_strc']['file']['name'])
        # try:
        #     f_cid_name = os_path.join(outdir, rrCache.__cache['cid_name']['file']['name'])
        # except KeyError:
        #     f_cid_name = ""
        #     logger.debug("   No cid_name file found in cache, skipping generation of cid_name")

        # Do not checksum since it is a dictionary
        if os_path.exists(f_cid_strc) and check_sha(
            f_cid_strc,
            rrCache.__cache['cid_strc']
        ): #and os_path.exists(f_cid_name) and check_sha(
        #     f_cid_name,
        #     rrCache.__cache['cid_name']
        # ):
            cid_strc = rrCache._load_json(f_cid_strc)
            logger.debug("   Cache file already exists")
        else:
            # if deprecatedCID_cid:
            #     if not deprecatedCID_cid['attr']:
            #         logger.debug("   Loading input data from file...")
            #         deprecatedCID_cid = rrCache._load_json(deprecatedCID_cid['file'])
            # else:
            #     deprecatedCID_cid = {'attr': {}}
            logger.debug("   Generating data...")
            # dep_files = [os_path.join(input_dir, 'metanetx', f) for f in rrCache.__cache['cid_strc']['deps']['file_deps']]
            dep_files = [os_path.join(input_dir, f) for f in rrCache.__cache['cid_strc']['deps']['file_deps']]
            dep_files += [os_path.join(input_dir, f) for f in rrCache.__cache['cid_xref']['deps']['file_deps']]
            cid_strc, cid_name = rrCache._m_mnxm_strc(dep_files, interactive=interactive, logger=logger)

            # if deprecatedCID_cid['attr'] != {}:
            #     # print(cid_strc['MNXM1106057'])
            #     # Replace compound IDs that have no structure with one that has.
            #     # Done from a manually built file
            #     with open(os_path.join(input_dir, 'MNXM_replacement_20190524.csv')) as csv_file:
            #         reader = csv_reader(csv_file, delimiter=' ')
            #         for row in reader:
            #             if not row[0].startswith('#') and len(row) > 1:
            #                 if row[1] != 'R_group':
            #                     # print(row)
            #                     cid_strc[row[0]] = cid_strc[row[1]]

            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(cid_strc, f_cid_strc, logger=logger)
            # if cid_name:
            #     rrCache._store_cache_to_file(cid_name, f_cid_name, logger=logger)
            # else:
            #     logger.debug("   No cid_name file found in cache, skipping generation of cid_name")

        return {
            'attr': cid_strc,
            'file': f_cid_strc
        }, {
            # 'attr': cid_name,
            # 'file': f_cid_name
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
                cid_strc['attr'] = rrCache._load_json(cid_strc['file'])
            logger.debug("   Generating data...")
            inchikey_cid = rrCache._m_inchikey_cid(cid_strc['attr'])
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(inchikey_cid, f_inchikey_cid, logger=logger)
            del inchikey_cid


    @staticmethod
    def _gen_cid_xref(
        input_dir: str,
        outdir: str,
        # deprecatedCID_cid: Dict,
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
            cid_xref = rrCache._load_json(f_cid_xref)
            logger.debug("   Cache file already exists")
        else:
            # if not deprecatedCID_cid['attr']:
            #     logger.debug("   Loading input data from file...")
            #     deprecatedCID_cid['attr'] = rrCache._load_json(deprecatedCID_cid['file'])
            logger.debug("   Generating data...")
            cid_xref = rrCache._m_mnxm_xref(
                os_path.join(input_dir, 'chem_xref.tsv'),
                # deprecatedCID_cid['attr']
            )
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(cid_xref, f_cid_xref, logger=logger)

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
            rrCache._store_cache_to_file(chebi_cid, f_chebi_cid, logger=logger)
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
            deprecatedRID_rid = rrCache._load_json(f_deprecatedRID_rid)
            logger.debug("   Cache file already exists")
        else:
            logger.debug("   Generating data...")
            deprecatedRID_rid = rrCache._m_deprecatedMNXR(
                os_path.join(input_dir, 'reac_xref.tsv')
            )
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(deprecatedRID_rid, f_deprecatedRID_rid, logger=logger)

        return {
            'attr': deprecatedRID_rid,
            'file': f_deprecatedRID_rid
        }


    @staticmethod
    def _gen_rr_reactions(
        input_dir: str,
        outdir: str,
        type: str = 'legacy',
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
            logger.debug("   Generating data...")
            # dep_files = [os_path.join(input_dir, 'metanetx', f) for f in rrCache.__cache[attribute]['deps']['file_deps']]
            dep_files = [os_path.join(input_dir, f) for f in rrCache.__cache[attribute]['deps']['file_deps']]
            if type == 'legacy':
                rr_reactions = rrCache._m_rr_reactions_legacy(
                    dep_files[0],
                    logger=logger
                )
            else:
                rr_reactions = rrCache._m_rr_reactions(
                    dep_files[0],
                    logger=logger
                )

            logger.debug("   Writing data to file...")

            rrCache._store_cache_to_file(rr_reactions, f_rr_reactions, logger=logger)

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
            rrCache._store_cache_to_file(comp_xref, f_comp_xref, logger=logger)
            # print_OK()
            del comp_xref
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(deprecatedCompID_compid, f_deprecatedCompID_compid, logger=logger)
            # print_OK()
            del deprecatedCompID_compid


    @staticmethod
    def _gen_template_reactions(
        input_dir: str,
        outdir: str,
        # deprecatedRID_rid: Dict,
        type: str = 'legacy',
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

        # dep_files = [os_path.join(input_dir, 'metanetx', f) for f in rrCache.__cache[attribute]['deps']['file_deps']]
        dep_files = [os_path.join(input_dir, f) for f in rrCache.__cache[attribute]['deps']['file_deps']]
        if type == 'legacy':
            template_reactions = rrCache._m_template_reactions_legacy(dep_files[0], logger=logger)
            # if deprecatedRID_rid:
            #     # Handle deprecated reaction IDs
            #     for depRID, newRID in deprecatedRID_rid['attr'].items():
            #         try:
            #             template_reactions[depRID] = template_reactions[newRID]
            #         except KeyError as key:
            #             logger.warning(f'Reaction ID {key} not found in {dep_files[0]}')
        else:
            template_reactions = rrCache._m_template_reactions(dep_files[0], logger=logger)

        logger.debug("   Writing data to file...")

        rrCache._store_cache_to_file(template_reactions, f_template_reactions, logger=logger)

        del template_reactions


    # def _load_from_file(self, filename):
    #     self.logger.debug("Loading "+filename+"...")
    #     data = self._load_json(
    #         os_path.join(self.__cache_dir, filename)
    #     )
    #     return data


    def _check_or_load_cache(self):
        self.logger.debug('Checking cache...')
        print_start(self.logger, 'Loading cache in memory')
        for attribute in self.__attributes_list:
            filename = rrCache.__cache[attribute]['file']['name']
            self.logger.debug('Loading '+attribute+' from '+filename+'...')
            if self.get(attribute) is None:
                self.set(
                    attribute,
                    self._load_json(
                        os_path.join(self.__cache_dir, filename)
                    )
                )
                dico = self._load_json(
                        os_path.join(self.__cache_dir, filename)
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
        logger.debug(f'Checking {file}...')

        # check if file exists
        if os_path.isfile(filename):
            logger.debug(f'File {file} present')
        else:
            logger.debug(f'File {filename} not present')
            # start_time = time_time()
            rrCache.__download_input_cache(url, file, outdir)
            # # if file size exceeds 20MB, compress it (compatible with all OS)
            # if os_path.getsize(filename) > 20 * 1024 * 1024:
            #     logger.debug(f'File {file} size exceeds 20MB, compressing...')
            #     with NamedTemporaryFile(delete=False) as temp_file:
            #         temp_file.close()
            #         # Compress the file
            #         with gzip_open(temp_file.name, 'wb') as f_out:
            #             with open(filename, 'rb') as f_in:
            #                 f_out.writelines(f_in)
            #     # Replace the original file with the compressed one
            #     os_rename(temp_file.name, filename+'.gz')
            #     logger.debug(f'File {file} compressed successfully')
            print_progress(logger)
            # end_time = time_time()

        # check sha512sum
        if check_sha(filename, fingerprint):
            logger.debug(f'File {file} present and sha ok')
        else:  # sha not ok
            logger.debug(f'\n\
                filename: {filename}\n\
                sha (computed): {sha512(Path(filename).read_bytes()).hexdigest()}\n\
                sha (expected): {fingerprint}\n\
            '
            )
            raise FileCorruptedError(f'Unable to download input-cache file {file}. Either the URL is broken or the file content has changed.')


    @staticmethod
    def __download_input_cache(
        url: str,
        file: str,
        outdir: str,
        logger: Logger = getLogger(__name__)
    ):

        """Download the input cache file from the given URL to the specified output directory.
        Args:
            url (str): The URL to download the file from.
            file (str): The name of the file to download.
            outdir (str): The directory to save the downloaded file.
        """
        logger.debug(f'Downloading {file} from {url}')

        # Check if the URL is empty
        if not url:
            logger.error(f'\n\n*** URL is empty, cannot download file {file}.\n')

        # Check if the output directory exists, if not create it
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
    def _load_json(filename, logger: Logger = getLogger(__name__)):
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
    def _store_cache_to_file(data, filename, logger: Logger = getLogger(__name__)):
        logger.debug(filename)
        if filename.endswith('.gz') or filename.endswith('.zip'):
            # Create the JSON string with sorted keys and no extra spaces
            # This ensures that the output is consistent for the same input data
            # which is important for reproducibility and caching
            json_bytes = json_dumps(data, sort_keys=True, separators=(",", ":")).encode("utf-8")
            # fp = gzip_open(filename, 'wt', encoding='ascii', mtime=0)
            # Write JSON into gzip file with reproducible output
            with open(filename, "wb") as raw:
                with GzipFile(fileobj=raw, mode="wb", mtime=0) as f:
                    f.write(json_bytes)
        else:
            try:
                fp = open(filename, 'w')
                json_dump(data, fp)
            except FileNotFoundError as e:
                logger.error(str(e))

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
        paths: List[str],
        # deprecatedCID_cid: Dict = None,
        interactive: bool = DEFAULTS['interactive'],
        logger: Logger = getLogger(__name__)
    ) -> Tuple[Dict, Dict]:

        """Parse the compounds.tsv file from RetroRules and the chem_prop.tsv file from MetanetX to generate a dictionary of compounds with their structures.
        Args:
            paths (List[str]): List of paths to the input files (rr_compounds_path first, then chem_prop_path).
            deprecatedCID_cid (Dict): Dictionary of deprecated CID to cid.
            interactive (bool): Whether to ask the user for confirmation before overwriting existing files.
        Returns:
            Tuple[Dict, Dict]: Tuple containing two dictionaries:
                - cid_strc: Dictionary of compounds with their structures.
                - cid_name: Dictionary of compound names.
        """

        logger.debug(f'paths: {paths}')
        # logger.debug(f'deprecatedCID_cid: {deprecatedCID_cid}')
        logger.debug(f'interactive: {interactive}')

        rr_compounds_path = paths[0]
        chem_prop_path = paths[1] if len(paths) > 1 else ""
        comp_xref_path = paths[2] if len(paths) > 2 else ""

        cid_strc = {}
        cid_name = {}

        # Parse the compounds.tsv file from RetroRules
        for row in csv_DictReader(gzip_open(rr_compounds_path, 'rt', encoding='utf-8-sig'), delimiter='\t'):
            if row.get('VALID', 'True').lower() != 'true':
                logger.debug('Skipping invalid compound entry: '+str(row))
                continue  # skip invalid entries
            row = {k.lower(): v for k, v in row.items()}  # normalize keys to lowercase
            if 'cid' not in row:
                # convert into 'id'
                row['cid'] = row.pop('id')
            tmp = {
                'formula':  row.get('formula', None),
                'inchi':    row.get('inchi', None),
                'inchikey': row.get('inchikey', None),
                'cid':      row.get('cid', row.get('ID', None)),
                'name':     row.get('name', None),
                'smiles':   row.get('smiles', row.get('SMILES', None))
            }
            logger.debug(f'Processing compound {tmp["cid"]} with InChI: {tmp["inchi"]} and InChIKey: {tmp["inchikey"]}')
            # for key in tmp.keys():
            #     try:
            #         tmp[key] = row[key]
            #     except KeyError:
            #         logger.debug(f'No {key} in RetroRules {rr_compounds_path} for '+str(row['cid'])+', setting to None')

            cid_strc[tmp['cid']] = tmp

        if chem_prop_path:
            # Parse the chem_prop.tsv file from MetanetX
            with open(chem_prop_path, 'rt') as f:
                tmp = {}
                # read CSV with both tab and space as delimiters
                c = csv_reader(f, delimiter='\t')
                for row in c:
                    if row[0].startswith('#'):
                        header = row
                    else:
                        # Set 'tmp' in a generic way
                        for i in range(len(header)):
                            # remove '#' from column field and
                            # convert to lower case
                            field = header[i].replace('#', '').lower()
                            if field == 'id':
                                field = 'cid'
                            tmp[field] = row[i]
                        if interactive:
                            print()
                            print("======================")
                            print()
                            print(tmp)
                        # mnxm = rrCache._checkCIDdeprecated(row[0], deprecatedCID_cid)
                        mnxm = row[0]
                        if interactive:
                            print(f'Converted into {mnxm}')
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
                        try: # in MetaNetX 3.1, there is no name field
                            if mnxm not in cid_name and tmp['name']:
                                cid_name[mnxm] = tmp['name']
                        except KeyError:
                            # If the name is not present
                            logger.debug('No name in chem_prop.tsv for '+str(mnxm)+', setting to None')
                        # Compound already in the dictionnary
                        if mnxm in cid_strc:
                            # # If the ID has been converted, then create a link
                            # if mnxm != row[0]:
                            #     cid_strc[mnxm] = cid_strc[row[0]]
                            if interactive:
                                print('already in cid_strc')
                            cid_strc[mnxm]['formula'] = row[2]
                            cid_strc[mnxm]['name'] = row[1]
                            for key in ['inchi', 'inchikey', 'smiles']:
                                if not cid_strc[mnxm][key] and tmp[key]:
                                    cid_strc[mnxm][key] = tmp[key]
                        # else:  # Compound not in the dictionnary
                        #     if interactive:
                        #         print('not yet in cid_strc')
                        #     # check to see if the inchikey is valid or not
                        #     otype = set({})
                        #     if not tmp['inchikey']:
                        #         otype.add('inchikey')
                        #     if not tmp['smiles']:
                        #         otype.add('smiles')
                        #     if not tmp['inchi']:
                        #         otype.add('inchi')
                        #     itype = ''
                        #     # Check if the current compound has
                        #     # InChI or SMILES description
                        #     if tmp['inchi']:
                        #         itype = 'inchi'
                        #     elif tmp['smiles']:
                        #         itype = 'smiles'
                        #     else:
                        #         ter = StreamHandler.terminator
                        #         StreamHandler.terminator = "\n"
                        #         logger.warning('No InChI or SMILES for '+str(tmp))
                        #         StreamHandler.terminator = ter
                        #     ter = StreamHandler.terminator
                        #     StreamHandler.terminator = "\n"
                        #     if itype:
                        #         try:
                        #             resConv = rrCache._convert_depiction(idepic=tmp[itype], itype=itype, otype=otype)
                        #             for i in resConv:
                        #                 tmp[i] = resConv[i]
                        #             logger.debug('Sructure conversion OK: '+str(tmp))
                        #         except rrCache.DepictionError as e:
                        #             logger.warning('Structure conversion FAILED: '+str(tmp))
                        #             logger.warning(e)
                        #     # StreamHandler.terminator = ter
                        #     # if mnxm == 'MNXM1106057':
                        #     #     print(tmp)
                        #     cid_strc[mnxm] = dict(tmp)
                        if interactive:
                            print(cid_strc[mnxm])
                            interactive = ask_user_input()
        else:
            logger.debug('No chem_prop.tsv or deprecatedCID_cid provided, skipping structure parsing')

        if comp_xref_path:
            # Add cross references from comp_xref.tsv
            comp_xref, _ = rrCache._m_mnxc_xref(
                comp_xref_path,
                logger=logger
            )
            for mnxc in comp_xref:
                if mnxc in cid_strc:
                    cid_strc[mnxc]['xref'] = comp_xref[mnxc]
                    logger.debug(f'Added cross references for {mnxc}: {comp_xref[mnxc]}')
        else:
            logger.debug('No comp_xref.tsv provided, skipping cross-reference parsing')
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
        # deprecatedCID_cid: Dict,
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        cid_xref = {}
        with open(chem_xref_path, 'rt') as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0] == '#':
                    # mnx = rrCache._checkCIDdeprecated(row[1], deprecatedCID_cid)
                    mnx = row[1]
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
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        rr_reactions = {}

        if not os_path.exists(rules_rall_path):
            logger.error('Could not read the rules file ('+str(rules_rall_path)+')')
            return None

        for row in csv_DictReader(gzip_open(rules_rall_path, 'rt'), delimiter='\t'):
            if row['TEMPLATE_ID'] not in rr_reactions:
                rr_reactions[row['TEMPLATE_ID']] = {}
            if row['REACTION_ID'] not in rr_reactions[row['TEMPLATE_ID']]:
                subtrates = {row['LEFT_IDS']: 1}
                products = dict(Counter(row['RIGHT_IDS'].split('.')))
                rr_reactions[row['TEMPLATE_ID']][row['REACTION_ID']] = {
                    'rule_id': row['TEMPLATE_ID'],
                    'rule_score': float(row['SCORE']),
                    'reac_id': row['REACTION_ID'],
                    'subs_id': row['LEFT_IDS'],
                    'rel_direction': (1 if row['DIRECTION'] == 'L2R' else -1),
                    'left': subtrates,
                    'right': products,
                    'left_excluded': row['LEFT_EXCLUDED_IDS'].split('.') if row['LEFT_EXCLUDED_IDS'] else [],
                    'right_excluded': row['RIGHT_EXCLUDED_IDS'].split('.') if row['RIGHT_EXCLUDED_IDS'] else []
                }
            # Handle multiple reactions per rule, update direction if needed
            else:
                if rr_reactions[row['TEMPLATE_ID']][row['REACTION_ID']]['rel_direction'] != 0 and (1 if row['DIRECTION'] == 'L2R' else -1) != rr_reactions[row['TEMPLATE_ID']][row['REACTION_ID']]['rel_direction']:
                    logger.debug('Updating direction for reaction '+str(row['REACTION_ID'])+' in rule '+str(row['TEMPLATE_ID'])+' from '+str(rr_reactions[row['TEMPLATE_ID']][row['REACTION_ID']]['rel_direction'])+' to bidirectional (0)')
                    rr_reactions[row['TEMPLATE_ID']][row['REACTION_ID']]['rel_direction'] = 0  # bidirectional

        return rr_reactions


    @staticmethod
    def _m_rr_reactions_legacy(
        rules_rall_path: str,
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        rr_reactions = {}

        if not os_path.exists(rules_rall_path):
            logger.error('Could not read the rules file ('+str(rules_rall_path)+')')
            return None

        for row in csv_DictReader(gzip_open(rules_rall_path, 'rt'), delimiter='\t'):
            # NOTE: as of now all the rules are generated using MNX
            # but it may be that other db are used, we are handling this case
            # WARNING: can have multiple products so need to seperate them
            products = {}

            for cid in row['Product_IDs'].split('.'):
                products[cid] = products.get(cid, 0) + 1

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
                    'reac_id': row['Reaction_ID'],
                    'subs_id': row['Substrate_ID'],
                    'rel_direction': int(row['Rule_relative_direction']),
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
        metadata_path: str,
        logger: Logger = getLogger(__name__)
    ) -> Dict:

        if not os_path.exists(metadata_path):
            logger.error('Cannot find file: '+str(metadata_path))
            return None

        reactions = {}

        for row in csv_DictReader(gzip_open(metadata_path, 'rt'), delimiter='\t'):
            if row['REACTION_ID'] not in reactions:
                # print(row)
                substrates = dict(Counter([row['LEFT_IDS']] + (row['LEFT_EXCLUDED_IDS'].split('.') if row['LEFT_EXCLUDED_IDS'] else [])))
                products = dict(Counter(row['RIGHT_IDS'].split('.') + (row['RIGHT_EXCLUDED_IDS'].split('.') if row['RIGHT_EXCLUDED_IDS'] else [])))
                # if row['REACTION_ID'] == 'MNXR182203':
                #     print(substrates)
                #     print(products)
                #     exit()
                main_left = row['LEFT_IDS']
                main_right = row['RIGHT_IDS'].split('.')[0]
                if row['DIRECTION'] == 'R2L':
                    # Swap left and right if direction is R2L
                    substrates, products = products, substrates
                    main_left, main_right = main_right, main_left
                    direction = -1
                else:
                    direction = 1
                reactions[row['REACTION_ID']] = {
                    'left': substrates,
                    'right': products,
                    'direction': direction,
                    'main_left': main_left,
                    'main_right': main_right
                }
            # Handle multiple reactions per rule, update direction if needed
            elif row['DIRECTION'] != reactions[row['REACTION_ID']]['direction']:
                reactions[row['REACTION_ID']]['direction'] = 0  # bidirectional

        return reactions


    @staticmethod
    def _m_template_reactions_legacy(
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
            if '@' in eq:
                regex_str = r'(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+'
            else:
                regex_str = r'(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)'
            for spe in re_findall(regex_str, eq.split('=')[side]):
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
    def _convert_depiction(idepic, itype='smiles', otype={'inchikey'}, logger=getLogger(__name__)):
        def MolFrom(idepic, itype, sanitize=True):
            if itype == 'smiles':
                return MolFromSmiles(idepic, sanitize=sanitize)
            elif itype == 'inchi':
                return MolFromInchi(idepic, sanitize=sanitize)
            else:
                raise NotImplementedError(
                    f'"{itype}" is not a valid input type for MolFrom'
                )
        # Import (if needed)
        rdmol = MolFrom(idepic, itype, sanitize=True)
        if rdmol is None:  # Check imprt
            logger.warning(
                f'Could not sanitize depiction "{idepic}" of type "{itype}", trying without sanitization'
            )
            rdmol = MolFrom(idepic, itype, sanitize=False)
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
