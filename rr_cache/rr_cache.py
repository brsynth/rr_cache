from os import (
    path as os_path,
    mkdir as os_mkdir,
    makedirs
)
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
from time       import time as time_time
from requests   import exceptions as r_exceptions
from redis      import StrictRedis
from hashlib    import sha512
from pathlib    import Path
from colored    import (
    attr as c_attr,
    fg as c_fg,
    bg as c_bg
)
from logging import (
    Logger,
    getLogger,
    StreamHandler
)
from brs_utils  import (
    print_OK,
    print_FAILED,
    print_start,
    print_progress,
    print_end,
    download,
)
from credisdict import (
    CRedisDict,
    wait_for_redis
)
from typing import (
    List,
    Tuple,
    Dict
)



#######################################################
################### rrCache  ##########################
#######################################################


class rrCache:
    """Class to generate the cache

    Contains all the functions that parse different files, used to calculate the thermodynamics and the FBA of the the other steps. These should be called only when the files have changes

    """

    # _input__cache_url = 'ftp://ftp.vital-it.ch/databases/metanetx/MNXref/3.2/'
    __cache_url       = 'https://gitlab.com/breakthewall/rpCache-data/-/raw/master/'

    # static attribues
    __convertMNXM = {
        'MNXM162231': 'MNXM6',
        'MNXM84':     'MNXM15',
        'MNXM96410':  'MNXM14',
        'MNXM114062': 'MNXM3',
        'MNXM145523': 'MNXM57',
        'MNXM57425':  'MNXM9',
        'MNXM137':    'MNXM588022'
        }

    # name: sha512sum
    __input__cache_files = {
            'chem_xref.tsv.gz':    'e558110990dcc75af943863790dc55360fd2d40ecb17d02335377671e80f0ab3738fd556acb340e03e48dd1afdec3eece1e92df1e18bc24e7445f24f778a10da',
            'MNXM_replacement_20190524.csv':    '0ffd5da832f9be057b19ca6a813a5ed3f53d7b6c247513f117450a63f5b34deb86de68ecbb5ed765ac6c3417f8976c4e39590a4dd18275351a076de863e6bfa9',
            'reac_xref.tsv.gz':    '48b991cf4a9c2ca573d395cf35c378881ed79e87772827647bfab2f6345499698664e07195ec10b342fc0164304dbd2363cccff1a1182225e6afebce3c16448b',
            'compounds.tsv.gz':    '719716bb880257bd014e045c03eb8dd12e2bbeba3aa52e38e9632ce605817b9dc09530e81fadd25542c0a439bdb81e1dfbd3a38f35b30b061845d1a880dbfe01',
            'chem_prop.tsv.gz':    'f2d220d1f0425e5e47f01e7deccfa46b60094d43b9f62b191ffb0fab8c00ef79e87c3b71d10bdcd26020608094f24884f51b3ebc3d7d3c9a6d594c6eaa324c66',
            'retrorules_rr02_flat_all.tsv.gz':   '890bdd24042c0192b5538964d775feefcb6cff9ad5f35690bfbfc5ae09334dd19df6828cdfc7f57a2018e090571517122b99d8760128052af898c638ae667e24',
            'comp_xref.tsv.gz':    '913a827f3645fda1699676ae6c32b9d7a8debae97ce7b0c386d8447f4eee5aa721d31bfb856d4092b3d5e987a8f19a6fe4bd28ddf1c5df5f85e71c3625bd1d81',
            'rxn_recipes.tsv.gz':  'dc0624f5ed7ab0b691d9a6ba02571a5cf334cfdb3109e78c98708e31574c46aeac2a97e9433788d80490ff80337679ccfd706cbb8e71a11cdc6122573bb69b0f'
            }

    # Attributes with dependencies (other attributes + input_cache files)
    __attributes_deps = {
            'deprecatedCID_cid': {
                'attr_deps': [],
                'file_deps': ['chem_xref.tsv.gz', 'MNXM_replacement_20190524.csv']
            },
            'deprecatedRID_rid': {
                'attr_deps': [],
                'file_deps': []
            },
            'cid_strc': {
                'attr_deps': ['deprecatedCID_cid'],
                'file_deps': ['compounds.tsv.gz', 'chem_prop.tsv.gz']
            },
            'cid_name': {
                'attr_deps': ['deprecatedCID_cid'],
                'file_deps': ['compounds.tsv.gz', 'chem_prop.tsv.gz']
            },
            'cid_xref': {
                'attr_deps': ['deprecatedCID_cid'],
                'file_deps': []
            },
            'chebi_cid': {
                'attr_deps': ['cid_xref'],
                'file_deps': []
            },
            'rr_reactions': {
                'attr_deps': ['deprecatedCID_cid', 'deprecatedRID_rid'],
                'file_deps': ['retrorules_rr02_flat_all.tsv.gz']
            },
            'inchikey_cid': {
                'attr_deps': ['cid_strc'],
                'file_deps': []
            },
            'comp_xref': {
                'attr_deps': [],
                'file_deps': ['comp_xref.tsv.gz']
            },
            'deprecatedCompID_compid': {
                'attr_deps': [],
                'file_deps': ['comp_xref.tsv.gz']
            },
            'rr_full_reactions': {
                'attr_deps': ['deprecatedCID_cid', 'deprecatedRID_rid'],
                'file_deps': ['rxn_recipes.tsv.gz']
            },
    }

    __attributes_list = list(__attributes_deps.keys())

    # name: sha512sum
    __cache_files = {
            __attributes_list[0]+'.json.gz': '1012d85b4720cf9706340e2d3bbef264dd95b67d510d8c1532c612a0f4aa6e1bdb2b77f36e135907070edb35038f31606f291ddd487df025e9999a4c543c739e',
            __attributes_list[1]+'.json.gz': '2cb1cc54c5a11962ef90a0da1d77e684b16b453214052541b473cf1b182491af9242d25e73f3d9b137e2c939eaf1ebf4adcee5b732ffab6cbb430d238f493b91',
            __attributes_list[2]+'.json.gz': '0139e8779c298fb089940ee3a57de0e4d13dcca40b954db1ca728dba0132fae43b7bf1adb128e6989433dae11c5c24a32f729aa7154c86e1e9236ac1b9dd5311',
            __attributes_list[3]+'.json.gz': 'fef4540c117f5ff75fe7106e55596b74729c8551d6713d454846a98509bc872114db92946389d57e24dcf6f588d85462e0d71606207c52676175fb5e5a81bf81',
            __attributes_list[4]+'.json.gz': '3c483178ea3d10a8b0ba9a5982cd70540af15fe4b95dd27844cad5c221554f253b75729dbd2ddd20b01dbe7f3b427e3adf13238c3fb1f3ae853feab9e4c36d9d',
            __attributes_list[5]+'.json.gz': '9579714acf9d25d6c207b94af5be9977a0d438910a5ef9c3cd22f10b3c3292097e2246dd7c3d387e61701d66d7fb25092f39639b12803cdd72494d04d5434158',
            __attributes_list[6]+'.json.gz': '70b20aa9cdf331c8407d86761596b851879d6878ecb7c8f89945c93d67569973845096c902c8b056064d16347fac380dfa17cd1918954a502d0fee70eca037aa',
            __attributes_list[7]+'.json.gz': 'd7d422e497af88bf8ea820857d4cca9846b2ccea7a9993beea479f270ba5852b901f3280a7c408d09ad8bd941cb2552940d281068c579592f41bbb659777a7b2',
            __attributes_list[8]+'.json.gz': '249a5bbd2b06c6326b3da8c759818c78bc13eed13b077c9f0ab1a2c912d01806e0e29d439fe3c09f5c5c6f87f7e02cd0bf0d7c5d01b81240aedb47b1bc5a664e',
            __attributes_list[9]+'.json.gz': '7e7e6a4805d74f680c31a48304abd4f221d394788604c935b91e977b78275cdf3266839d8bc1097637c2fd633b3583f082b60f0d5a54985bba931f0f945c9802',
            __attributes_list[10]+'.json.gz': '33e6f15a61352686169d86c8c444717eb6ea315b4683ca1d60a8b0d542a1fea2411fc9327d2604df1e47defa836616e60beb1b139a43249cb44bb9b0e91a1c70'
            }
    # __cache_files = {
    #         __attributes_list[0]+'.json.gz': '698a3e83cf4f9206ea2644c9c35a9af53957838baaae6efb245d02b6b8d0ea8b25c75008e562b99ba3e0189e50ee47655376f2d0635f6206e0015f91f0e4bad8',
    #         __attributes_list[1]+'.json.gz': '51554c6f6ae99c6755da7496208b3feec30547bc4cf3007d9fd30f46fa4c0cc73bad5aeb743dca07e32711c4346504296bee776d135fb18e96c891a0086fc87e',
    #         __attributes_list[2]+'.json.gz': '0021ef63165d75ee6b8c209ccf14b8a1b8b7b263b4077f544729c47b5525f66511c3fa578fd2089201abb61693085b9912639e62f7b7481d06ad1f38bfc2dd8e',
    #         __attributes_list[3]+'.json.gz': '7d559cc7389c0cb2bd10f92e6e845bb5724be64d1624adc4e447111fc63599bb69396cd0cc3066a6bb19910c00e266c97e21b1254d9a6dc9da3a8b033603fcff',
    #         __attributes_list[4]+'.json.gz': '587d6c5206ee94e63af6d9eaf49fd5e2ca417308b3ece8a7f47e916c42376e2c8635a031ce26dc815cd7330f2323054a44d23951e416a9a29c5a9a2ab51e8953',
    #         __attributes_list[5]+'.json.gz': '8783aaa65a281c4a7ab3a82a6dc99620418ed2be4a739f46db8ee304fcb3536a78fed5a955e1c373a20c3e7d3673793157c792b4429ecb5c68ddaddb1a0f7de7',
    #         __attributes_list[6]+'.json.gz': '8007480fc607caf41f0f9a93beb66c7caa66c37a3d01a809f6b94bc0df469cec72091e8cc0fbabb3bd8775e9776b928ecda2779fc545c7e4b9e71c504f9510ce',
    #         __attributes_list[7]+'.json.gz': 'afc2ad3d31366a8f7fe1604fa49c190ade6d46bc8915f30bd20fdfdfc663c979bb10ca55ad10cadec6002a17add46639c70e7adf89cb66c57ed004fd3e4f0051',
    #         __attributes_list[8]+'.json.gz': '81c673fe1940e25a6a9722fd74b16bc30e1590db0c40810f541ad4ffba7ae04c01268b929d4bf944e84095a0c2a1d0079d1861bc1df3e8308fbb6b35e0aaf107',
    #         __attributes_list[9]+'.json.gz': '599e4de4935d2ba649c0b526d8aeef6f0e3bf0ed9ee20adad65cb86b078ac139e4cc9758945c2bb6da1c6840867239c5415cb5bceeb80164798ff627aac0a985',
    #         __attributes_list[10]+'.json.gz': '599e4de4935d2ba649c0b526d8aeef6f0e3bf0ed9ee20adad65cb86b078ac139e4cc9758945c2bb6da1c6840867239c5415cb5bceeb80164798ff627aac0a985'
    #         }

    __ext = '.json.gz'


    ## Cache constructor
    #
    # @param self The object pointer
    # @param db Mode of storing objects ('file' or 'redis')
    def __init__(self, db='file', attrs=[], logger=getLogger(__name__)):

        self.logger = logger

        self.logger.debug('New instance of rrCache')

        self.store_mode = db
        rrCache._db_timeout = 10

        self.dirname = os_path.dirname(os_path.abspath( __file__ ))#+"/.."
        # # input_cache
        # self._input__cache_dir = os_path.join(self.dirname, 'input_cache')
        # cache
        self.__cache_dir = os_path.join(self.dirname, 'cache')

        # if self.__attributes_list is not None:
        self.load(attrs)


    def load(self, attrs: List=[]):

        if attrs is None:
            return

        if attrs != []:
            if not isinstance(attrs, list):
                self.logger.warning('\'attrs\' argument is not of type list, trying to convert...')
                self.__attributes_list = [attrs]
                self.logger.warning('OK')
            else:
                self.__attributes_list = attrs

        if self.store_mode!='file':
            self.redis = StrictRedis(host=self.store_mode, port=6379, db=0, decode_responses=True)
            if not wait_for_redis(self.redis, self._db_timeout):
                self.logger.critical("Database "+self.store_mode+" is not reachable")
                exit()
            for attr in self.__attributes_list:
                setattr(self, '__'+attr, CRedisDict(attr, self.redis))
        else:
            for attr in self.__attributes_list:
                setattr(self, '__'+attr, None)

        rrCache._check_or_download_cache_to_disk(
            self.__cache_dir,
            self.__attributes_list,
            self.logger
        )
        # rrCache.generate_cache(self.__cache_dir)
        try:
            self._check_or_load_cache()
        except (r_exceptions.RequestException,
                r_exceptions.InvalidSchema,
                r_exceptions.ConnectionError):
            rrCache.generate_cache(self.__cache_dir)
            self._check_or_load_cache()
        # exit()
        # try:
        #     self._check_or_load_cache()
        # except FileNotFoundError:
        #     try:
        #         StreamHandler.terminator = "\r"
        #         logger.info('')
        #         rrCache._check_or_download_cache_to_disk(
        #             self.__cache_dir,
        #             self.__attributes_list,
        #             self.logger
        #         )
        #         self._check_or_load_cache()
        #     except (r_exceptions.RequestException,
        #             r_exceptions.InvalidSchema,
        #             r_exceptions.ConnectionError):
        #         print_FAILED()
        #         rrCache.generate_cache(self.__cache_dir)
        #         self._check_or_load_cache()


    @staticmethod
    def _check_or_download_cache_to_disk(
        cache_dir: str,
        attributes: List,
        logger: Logger=getLogger(__name__)
    ):
        logger.debug('cache_dir: '+str(cache_dir))
        logger.debug('attributes: '+str(attributes))

        print_start(logger, 'Downloading cache')

        for attr in attributes:
            print_progress()
            filename = attr+rrCache.__ext
            if os_path.isfile(
                os_path.join(cache_dir, filename)
            ) and sha512(
                Path(
                    os_path.join(cache_dir, filename)
                ).read_bytes()
            ).hexdigest() == rrCache.__cache_files[filename]:
                logger.debug(filename+" already downloaded")
                # print_OK()
            else:
                filename = attr+rrCache.__ext
                logger.debug("Downloading "+filename+"...")
                start_time = time_time()
                if not os_path.isdir(cache_dir):
                    os_mkdir(cache_dir)
                download(
                    rrCache.__cache_url+filename,
                    os_path.join(cache_dir, filename)
                )
                rrCache.__cache_files[attr] = True
                end_time = time_time()
                # print_OK(end_time-start_time)

        print_end(logger)


    def get(self, attr: str):
        try:
            return getattr(self, '__'+attr)
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
            #self.expression = expression
            self.message = message

    #url = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/'
    #url = 'ftp://ftp.vital-it.ch/databases/metanetx/MNXref/3.2/'

    @staticmethod
    def generate_cache(outdir, logger=getLogger(__name__)):

        if outdir == '':
            outdir = 'cache'
            input_dir = 'input-'+outdir
        else:
            input_dir = os_path.join(
                outdir,
                'input-cache'
            )
            outdir = os_path.join(
                outdir,
                'cache'
            )
        if not os_path.isdir(outdir):
            makedirs(outdir)
        if not os_path.isdir(input_dir):
            makedirs(input_dir)

        # FETCH INPUT_CACHE FILES
        url = rrCache.__cache_url
        # input_dir = os_path.join(
        #     os_path.normpath(outdir),
        # )'input-'+os_path.basename(os_path.normpath(outdir))
        print_start(logger, 'Downloading input cache')
        for file in rrCache.__input__cache_files.keys():
            rrCache._download_input_cache(url, file, input_dir)
            print_progress(logger)
        print_end(logger)

        # GENERATE CACHE FILES AND STORE THEM TO DISK
        print_start(logger, 'Generating cache')
        deprecatedCID_cid  = rrCache._gen_deprecatedCID_cid(input_dir, outdir, logger)
        print_progress(logger)
        cid_strc, cid_name = rrCache._gen_cid_strc_cid_name(input_dir, outdir, deprecatedCID_cid, logger)
        print_progress(logger)
        rrCache._gen_inchikey_cid(input_dir, outdir, cid_strc, logger)
        print_progress(logger)
        del cid_strc, cid_name
        cid_xref           = rrCache._gen_cid_xref(input_dir, outdir, deprecatedCID_cid, logger)
        print_progress(logger)
        rrCache._gen_chebi_cid(input_dir, outdir, cid_xref)
        print_progress(logger)
        del cid_xref
        deprecatedRID_rid  = rrCache._gen_deprecatedRID_rid(input_dir, outdir, logger)
        print_progress(logger)
        rrCache._gen_rr_reactions(input_dir, outdir, deprecatedCID_cid, deprecatedRID_rid, logger)
        print_progress(logger)
        rrCache._gen_comp_xref_deprecatedCompID_compid(input_dir, outdir, logger)
        print_progress(logger)
        rrCache._gen_rr_full_reactions(input_dir, outdir, deprecatedCID_cid, deprecatedRID_rid, logger)
        print_progress(logger)
        del deprecatedCID_cid, deprecatedRID_rid
        print_progress(logger)
        print_end(logger)


    @staticmethod
    def _gen_deprecatedCID_cid(
        input_dir: str,
        outdir: str,
        logger=getLogger(__name__)
    ) -> Dict:
        attribute = 'deprecatedCID_cid'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        deprecatedCID_cid = None
        f_deprecatedCID_cid = os_path.join(outdir, attribute)+rrCache.__ext

        if not os_path.isfile(f_deprecatedCID_cid):
            logger.debug("   Generating data...")
            deprecatedCID_cid = rrCache._m_deprecatedMNXM(
                os_path.join(input_dir, 'chem_xref.tsv.gz')
            )
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(deprecatedCID_cid, f_deprecatedCID_cid)

        else:
            deprecatedCID_cid = rrCache._load_cache_from_file(f_deprecatedCID_cid)
            logger.debug("   Cache file already exists")

        return {
            'attr': deprecatedCID_cid,
            'file': f_deprecatedCID_cid
        }


    @staticmethod
    def _gen_cid_strc_cid_name(
        input_dir: str,
        outdir: str,
        deprecatedCID_cid: Dict,
        logger=getLogger(__name__)
    ) -> Dict:

        attribute = 'cid_strc, cid_name'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        cid_strc = None
        cid_name = None
        f_cid_strc = os_path.join(outdir, 'cid_strc')+rrCache.__ext
        f_cid_name = os_path.join(outdir, 'cid_name')+rrCache.__ext

        if not os_path.isfile(f_cid_strc):
            if not deprecatedCID_cid['attr']:
                logger.debug("   Loading input data from file...")
                deprecatedCID_cid = rrCache._load_cache_from_file(deprecatedCID_cid['file'])
                print_OK()
            logger.debug("   Generating data...")
            cid_strc, cid_name = rrCache._m_mnxm_strc(
                os_path.join(input_dir, 'compounds.tsv.gz'),
                os_path.join(input_dir, 'chem_prop.tsv.gz'),
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

        else:
            cid_strc = rrCache._load_cache_from_file(f_cid_strc)
            logger.debug("   Cache file already exists")
            
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
        logger=getLogger(__name__)
    ) -> None:
        attribute = 'inchikey_cid'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        inchikey_cid = None
        f_inchikey_cid = os_path.join(outdir, attribute)+rrCache.__ext
        if not os_path.isfile(f_inchikey_cid):
            if not cid_strc['attr']:
                logger.debug("   Loading input data from file...")
                cid_strc['attr'] = rrCache._load_cache_from_file(cid_strc['file'])
            logger.debug("   Generating data...")
            inchikey_cid = rrCache._m_inchikey_cid(cid_strc['attr'])
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(inchikey_cid, f_inchikey_cid)
        else:
            logger.debug("   Cache file already exists")


    @staticmethod
    def _gen_cid_xref(
        input_dir: str,
        outdir: str,
        deprecatedCID_cid: Dict,
        logger=getLogger(__name__)
    ) -> Dict:
        attribute = 'cid_xref'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        cid_xref = None
        f_cid_xref = os_path.join(outdir, attribute)+rrCache.__ext

        if not os_path.isfile(f_cid_xref):
            if not deprecatedCID_cid['attr']:
                logger.debug("   Loading input data from file...")
                deprecatedCID_cid['attr'] = rrCache._load_cache_from_file(deprecatedCID_cid['file'])
            logger.debug("   Generating data...")
            cid_xref = rrCache._m_mnxm_xref(
                os_path.join(input_dir, 'chem_xref.tsv.gz'),
                deprecatedCID_cid['attr']
            )
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(cid_xref, f_cid_xref)

        else:
            cid_xref = rrCache._load_cache_from_file(f_cid_xref)
            logger.debug("   Cache file already exists")

        return {
            'attr': cid_xref,
            'file': f_cid_xref
        }


    @staticmethod
    def _gen_chebi_cid(
        input_dir: str,
        outdir: str,
        cid_xref: Dict,
        logger=getLogger(__name__)
    ) -> Dict:
        attribute = 'chebi_cid'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        chebi_cid = None
        f_chebi_cid = os_path.join(outdir, attribute)+rrCache.__ext
        if not os_path.isfile(f_chebi_cid):
            logger.debug("   Generating data...")
            chebi_cid = rrCache._m_chebi_cid(cid_xref['attr'])
            # print_OK()
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(chebi_cid, f_chebi_cid)
            del chebi_cid
            # print_OK()
        else:
            logger.debug("   Cache file already exists")
            # print_OK()


    @staticmethod
    def _gen_deprecatedRID_rid(
        input_dir: str,
        outdir: str,
        logger=getLogger(__name__)
    ) -> Dict:
        attribute = 'deprecatedRID_rid'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        deprecatedRID_rid = None
        f_deprecatedRID_rid = os_path.join(outdir, attribute)+rrCache.__ext
        if not os_path.isfile(f_deprecatedRID_rid):
            logger.debug("   Generating data...")
            deprecatedRID_rid = rrCache._m_deprecatedMNXR(
                os_path.join(input_dir, 'reac_xref.tsv.gz')
            )
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(deprecatedRID_rid, f_deprecatedRID_rid)
        else:
            deprecatedRID_rid = rrCache._load_cache_from_file(f_deprecatedRID_rid)
            logger.debug("   Cache file already exists")
        return {
            'attr': deprecatedRID_rid,
            'file': f_deprecatedRID_rid
        }


    @staticmethod
    def _gen_rr_reactions(
        input_dir: str,
        outdir: str,
        deprecatedCID_cid: Dict,
        deprecatedRID_rid: Dict,
        logger=getLogger(__name__)
    ) -> None:
        attribute = 'rr_reactions'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        rr_reactions = None
        f_rr_reactions = os_path.join(outdir, attribute)+rrCache.__ext
        if not os_path.isfile(f_rr_reactions):
            if not deprecatedCID_cid['attr']:
                logger.debug("   Loading input data from file...")
                deprecatedCID_cid['attr'] = rrCache._load_cache_from_file(deprecatedCID_cid['file'])
            if not deprecatedRID_rid['attr']:
                logger.debug("   Loading input data from file...")
                deprecatedRID_rid['attr'] = rrCache._load_cache_from_file(deprecatedRID_rid['file'])
                print_OK()
            logger.debug("   Generating data...")
            rr_reactions = rrCache._m_rr_reactions(
                os_path.join(input_dir, 'retrorules_rr02_flat_all.tsv.gz'),
                deprecatedCID_cid,
                deprecatedRID_rid
            )
            del deprecatedRID_rid
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(rr_reactions, f_rr_reactions)
            del rr_reactions
        else:
            logger.debug("   Cache file already exists")


    @staticmethod
    def _gen_comp_xref_deprecatedCompID_compid(
        input_dir: str,
        outdir: str,
        logger=getLogger(__name__)
    ) -> None:
        attribute = 'comp_xref, deprecatedCompID_compid'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        comp_xref = deprecatedCompID_compid = None
        f_comp_xref = os_path.join(outdir, 'comp_xref')+rrCache.__ext
        f_deprecatedCompID_compid = outdir+'deprecatedCompID_compid'+rrCache.__ext
        if not os_path.isfile(f_comp_xref) or not os_path.isfile(f_deprecatedCompID_compid):
            logger.debug("   Generating data...")
            comp_xref,deprecatedCompID_compid = rrCache._m_mnxc_xref(
                os_path.join(input_dir, 'comp_xref.tsv.gz')
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
        else:
            logger.debug("   Cache files already exist")
            # print_OK()


    @staticmethod
    def _gen_rr_full_reactions(
        input_dir: str,
        outdir: str,
        deprecatedCID_cid: Dict,
        deprecatedRID_rid: Dict,
        logger=getLogger(__name__)
    ) -> None:
        attribute = 'rr_full_reactions'
        logger.debug(c_attr('bold')+attribute+c_attr('reset'))
        rr_full_reactions = None
        f_rr_full_reactions = os_path.join(outdir, attribute)+rrCache.__ext
        if not os_path.isfile(f_rr_full_reactions):
            logger.debug("   Generating data...")
            if not deprecatedCID_cid['attr']:
                logger.debug("   Loading input data from file...")
                deprecatedCID_cid = rrCache._load_cache_from_file(deprecatedCID_cid['file'])
            if not deprecatedRID_rid:
                logger.debug("   Loading input data from file...")
                deprecatedRID_rid = rrCache._load_cache_from_file(deprecatedRID_rid['file'])
            rr_full_reactions = rrCache._m_rr_full_reactions(
                os_path.join(input_dir, 'rxn_recipes.tsv.gz'),
                deprecatedCID_cid['attr'],
                deprecatedRID_rid['attr']
            )
            logger.debug("   Writing data to file...")
            rrCache._store_cache_to_file(rr_full_reactions, f_rr_full_reactions)
            del rr_full_reactions
        else:
            logger.debug("   Cache file already exists")


    def _load_from_file(self, attribute):
        filename = attribute+rrCache.__ext
        self.logger.debug("Loading "+filename+"...")
        data = self._load_cache_from_file(
            os_path.join(self.__cache_dir, filename)
        )
        return data


    def _check_or_load_cache(self):
        self.logger.debug('store_mode: '+self.store_mode)
        if self.store_mode=='file':
            self._check_or_load_cache_in_memory()
        else:
            self._check_or_load_cache_in_db()


    def _check_or_load_cache_in_memory(self):
        print_start(self.logger, 'Loading cache in memory')
        for attribute in self.__attributes_list:
            if self.get(attribute) is None:
                self.set(
                    attribute,
                    self._load_from_file(attribute)
                )
                print_progress(self.logger)
            else:
                self.logger.debug(attribute+" already loaded in memory")
        print_end(self.logger)

    def _check_or_load_cache_in_db(self):
        print_start(self.logger, 'Loading cache in db')
        for attribute in self.__attributes_list:
            if not CRedisDict.exists(self.redis, attribute):
                self._store_cache_to_db(attribute, self._load_from_file(attribute))
            else:
                self.logger.debug(attribute+" already loaded in db")
            print_progress(self.logger)
        print_end(self.logger)


    @staticmethod
    def _download_input_cache(url, file, outdir):
        if not os_path.isdir(outdir):
            os_mkdir(outdir)
        filename = os_path.join(outdir, file)
        if not os_path.isfile(filename):
            start_time = time_time()
            rrCache.__download_input_cache(url, file, outdir)
            end_time = time_time()


    @staticmethod
    def __download_input_cache(url, file, outdir):

        if not os_path.isdir(outdir):
            os_mkdir(outdir)

        if file in [
            'reac_xref.tsv.gz',
            'chem_xref.tsv.gz',
            'chem_prop.tsv.gz',
            'comp_xref.tsv.gz'
        ]:
            download(
                url+'metanetx/'+file,
                os_path.join(outdir, file)
            )

        elif file in [
            'compounds.tsv.gz',
            'rxn_recipes.tsv.gz'
        ]:
            download(
                url+'rr02_more_data/'+file,
                os_path.join(outdir, file)
            )

        elif file == 'retrorules_rr02_flat_all.tsv.gz':
            download(
                url+'retrorules_rr02_rp3_hs/'+file,
                os_path.join(outdir, file)
            )

        else:
            download(
                url+file,
                os_path.join(outdir, file)
            )


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
    def _load_cache_from_file(filename):
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

    ## Method to store data into redis database
    #
    #  Assign a CRedisDict object to the attribute to copy data into the database
    #
    #  @param self Object pointer
    #  @param attr_name Attribute name (database key)
    #  @param data Content of the attribute
    def _store_cache_to_db(self, attr_name, data):
        print("Storing "+attr_name+" to db...", end = '', flush=True)
        setattr(rrCache, attr_name, CRedisDict(attr_name, self.redis, data))
        print_OK()

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
        with gzip_open(xref_path, 'rt') as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = row[0].split(':')
                    if mnx[0]=='deprecated':
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
    def _m_mnxm_strc(rr_compounds_path, chem_prop_path, deprecatedCID_cid, logger=getLogger(__name__)):

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
                resConv = rrCache._convert_depiction(idepic=tmp['inchi'], itype='inchi', otype={'smiles','inchikey'})
                for i in resConv:
                    tmp[i] = resConv[i]
            except rrCache.DepictionError as e:
                logger.warning('Could not convert some of the structures: '+str(tmp))
                logger.warning(e)
            cid_strc[tmp['cid']] = tmp
        with gzip_open(chem_prop_path, 'rt') as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnxm = rrCache._checkCIDdeprecated(row[0], deprecatedCID_cid)
                    tmp = {
                        'formula':  row[2],
                        'smiles': row[6],
                        'inchi': row[5],
                        'inchikey': row[8],
                        'cid': mnxm,
                        'name': row[1]
                    }
                    for i in tmp:
                        if tmp[i]=='' or tmp[i]=='NA':
                            tmp[i] = None
                    if not mnxm in cid_name and tmp['name']:
                        cid_name[mnxm] = tmp['name']
                    if mnxm in cid_strc:
                        cid_strc[mnxm]['formula'] = row[2]
                        cid_strc[mnxm]['name'] = row[1]
                        if not cid_strc[mnxm]['smiles'] and tmp['smiles']:
                            cid_strc[mnxm]['smiles'] = tmp['smiles']
                        if not cid_strc[mnxm]['inchikey'] and tmp['inchikey']:
                            cid_strc[mnxm]['inchikey'] = tmp['inchikey']
                    else:
                        #check to see if the inchikey is valid or not
                        otype = set({})
                        if not tmp['inchikey']:
                            otype.add('inchikey')
                        if not tmp['smiles']:
                            otype.add('smiles')
                        if not tmp['inchi']:
                            otype.add('inchi')
                        itype = ''
                        if tmp['inchi']:
                            itype = 'inchi'
                        elif tmp['smiles']:
                            itype = 'smiles'
                        else:
                            ter = StreamHandler.terminator
                            StreamHandler.terminator = "\n"
                            logger.warning('No valid entry for the convert_depiction function')
                            StreamHandler.terminator = ter
                            continue
                        try:
                            resConv = rrCache._convert_depiction(idepic=tmp[itype], itype=itype, otype=otype)
                            for i in resConv:
                                tmp[i] = resConv[i]
                        except rrCache.DepictionError as e:
                            ter = StreamHandler.terminator
                            StreamHandler.terminator = "\n"
                            logger.warning('Could not convert some of the structures: '+str(tmp))
                            logger.warning(e)
                            StreamHandler.terminator = ter
                        cid_strc[tmp['cid']] = tmp
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
    def _m_mnxm_xref(chem_xref_path, deprecatedCID_cid, logger=getLogger(__name__)):
        cid_xref = {}
        with gzip_open(chem_xref_path, 'rt') as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = rrCache._checkCIDdeprecated(row[1], deprecatedCID_cid)
                    if len(row[0].split(':'))==1:
                        dbName = 'mnx'
                        dbId = row[0]
                    else:
                        dbName = row[0].split(':')[0]
                        dbId = ''.join(row[0].split(':')[1:])
                        if dbName=='deprecated':
                            dbName = 'mnx'
                    #mnx
                    if not mnx in cid_xref:
                        cid_xref[mnx] = {}
                    if not dbName in cid_xref[mnx]:
                        cid_xref[mnx][dbName] = []
                    if not dbId in cid_xref[mnx][dbName]:
                        cid_xref[mnx][dbName].append(dbId)
                    ### DB ###
                    if not dbName in cid_xref:
                        cid_xref[dbName] = {}
                    if not dbId in cid_xref[dbName]:
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
    def _m_mnxc_xref(comp_xref_path, logger=getLogger(__name__)):
        comp_xref = {}
        deprecatedCompID_compid = {}
        try:
            with gzip_open(comp_xref_path, 'rt') as f:
                c = csv_reader(f, delimiter='\t')
                #not_recognised = []
                for row in c:
                    #cid = row[0].split(':')
                    if not row[0][0]=='#':
                        #collect the info
                        mnxc = row[1]
                        if len(row[0].split(':'))==1:
                            dbName = 'mnx'
                            dbCompId = row[0]
                        else:
                            dbName = row[0].split(':')[0]
                            dbCompId = ''.join(row[0].split(':')[1:])
                            dbCompId = dbCompId.lower()
                        if dbName=='deprecated':
                            dbName = 'mnx'
                        #create the dicts
                        if not mnxc in comp_xref:
                            comp_xref[mnxc] = {}
                        if not dbName in comp_xref[mnxc]:
                            comp_xref[mnxc][dbName] = []
                        if not dbCompId in comp_xref[mnxc][dbName]:
                            comp_xref[mnxc][dbName].append(dbCompId)
                        #create the reverse dict
                        if not dbCompId in deprecatedCompID_compid:
                            deprecatedCompID_compid[dbCompId] = mnxc
        except FileNotFoundError:
            logger.error('comp_xref file not found')
            return {}
        return comp_xref,deprecatedCompID_compid


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
    def _m_rr_reactions(rules_rall_path, deprecatedCID_cid, deprecatedRID_rid, logger=getLogger(__name__)):
        rr_reactions = {}
        try:
            for row in csv_DictReader(gzip_open(rules_rall_path, 'rt'), delimiter='\t'):
                # NOTE: as of now all the rules are generated using MNX
                # but it may be that other db are used, we are handling this case
                # WARNING: can have multiple products so need to seperate them
                products = {}
                for i in row['Product_IDs'].split('.'):
                    cid = rrCache._checkCIDdeprecated(i, deprecatedCID_cid)
                    if not cid in products:
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
                        'reac_id': rrCache._checkRIDdeprecated(row['Reaction_ID'], deprecatedRID_rid),
                        'subs_id': rrCache._checkCIDdeprecated(row['Substrate_ID'], deprecatedCID_cid),
                        'rel_direction': int(row['Rule_relative_direction']),
                        'left': {rrCache._checkCIDdeprecated(row['Substrate_ID'], deprecatedCID_cid): 1},
                        'right': products}
                except ValueError:
                    logger.error('Problem converting rel_direction: '+str(row['Rule_relative_direction']))
                    logger.error('Problem converting rule_score: '+str(row['Score_normalized']))
            return rr_reactions
        except FileNotFoundError as e:
                logger.error('Could not read the rules_rall file ('+str(rules_rall_path)+')')
                return {}


    ## Generate complete reactions from the rxn_recipes.tsv from RetroRules
    #
    #  These are the compplete reactions from which the reaction rules are generated from. This is used to
    #  reconstruct the full reactions from monocomponent reactions
    #  Structur of the return: rr_full_reactions['MNXR142257'] = {'left': {'MNXM4660': 1}, 'right': {'MNXM97172': 1}, 'direction': 0, 'main_left': ['MNXM4660'], 'main_right': ['MNXM97172']}
    #
    #  @param self The pointer object
    #  @param rxn_recipes_path Path to the recipes file
    #  @return Boolean that determines the success or failure of the function
    @staticmethod
    def _m_rr_full_reactions(rxn_recipes_path, deprecatedCID_cid, deprecatedRID_rid, logger=getLogger(__name__)):
        #### for character matching that are returned
        DEFAULT_STOICHIO_RESCUE = {"4n": 4, "3n": 3, "2n": 2, 'n': 1,
                           '(n)': 1, '(N)': 1, '(2n)': 2, '(x)': 1,
                           'N': 1, 'm': 1, 'q': 1,
                           '0.01': 1, '0.1': 1, '0.5': 1, '1.5': 1,
                           '0.02': 1, '0.2': 1,
                           '(n-1)': 0, '(n-2)': -1}
        reaction = {}

        try:

            for row in csv_DictReader(gzip_open(rxn_recipes_path, 'rt'), delimiter='\t'):
                tmp = {} # makes sure that if theres an error its not added
                # parse the reaction equation
                if not len(row['Equation'].split('='))==2:
                    logger.warning('There should never be more or less than a left and right of an equation')
                    logger.warnin(row['Equation'])
                    continue

                ######### LEFT ######
                #### MNX id
                tmp['left'] = {}
                for spe in re_findall(r'(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+', row['Equation'].split('=')[0]):
                    # 1) try to rescue if its one of the values
                    try:
                        tmp['left'][rrCache._checkCIDdeprecated(spe[1], deprecatedCID_cid)] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                    except KeyError:
                        # 2) try to convert to int if its not
                        try:
                            tmp['left'][rrCache._checkCIDdeprecated(spe[1], deprecatedCID_cid)] = int(spe[0])
                        except ValueError:
                            ter = StreamHandler.terminator
                            StreamHandler.terminator = "\n"
                            logger.warning('Cannot convert '+str(spe[0]))
                            StreamHandler.terminator = ter
                            continue

                ####### RIGHT #####
                #### MNX id
                tmp['right'] = {}
                for spe in re_findall(r'(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+', row['Equation'].split('=')[1]):
                    # 1) try to rescue if its one of the values
                    try:
                        tmp['right'][rrCache._checkCIDdeprecated(spe[1], deprecatedCID_cid)] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                    except KeyError:
                        # 2) try to convert to int if its not
                        try:
                            tmp['right'][rrCache._checkCIDdeprecated(spe[1], deprecatedCID_cid)] = int(spe[0])
                        except ValueError:
                            ter = StreamHandler.terminator
                            StreamHandler.terminator = "\n"
                            logger.warning('Cannot convert '+str(spe[0]))
                            StreamHandler.terminator = ter
                            continue

                ####### DIRECTION ######
                try:
                    tmp['direction'] = int(row['Direction'])
                except ValueError:
                    logger.error('Cannot convert '+str(row['Direction'])+' to int')
                    continue
                ### add the others
                tmp['main_left'] = row['Main_left'].split(',')
                tmp['main_right'] = row['Main_right'].split(',')

                rid = row['#Reaction_ID']
                new_rid = rrCache._checkRIDdeprecated(row['#Reaction_ID'], deprecatedRID_rid)
                reaction[new_rid] = tmp
                if new_rid != rid: # rid is deprecated and has to point toward new_rid
                    reaction[rid] = reaction[new_rid]
 
            return reaction

        except FileNotFoundError:
            logger.error('Cannot find file: '+str(rxn_recipes_path))
            return {}


    ######################## Generic functions ###############################

    ## Convert chemical depiction to others type of depictions
    #
    # Usage example:
    # - convert_depiction(idepic='CCO', otype={'inchi', 'smiles', 'inchikey'})
    # - convert_depiction(idepic='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', itype='inchi', otype={'inchi', 'smiles', 'inchikey'})
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
            raise NotImplementedError('"{}" is not a valid input type'.format(itype))
        if rdmol is None:  # Check imprt
            raise rrCache.DepictionError('Import error from depiction "{}" of type "{}"'.format(idepic, itype))
        # Export
        odepic = dict()
        for item in otype:
            if item == 'smiles':
                odepic[item] = MolToSmiles(rdmol)  # MolToSmiles is tricky, one mays want to check the possible options..
            elif item == 'inchi':
                odepic[item] = MolToInchi(rdmol)
            elif item == 'inchikey':
                odepic[item] = MolToInchiKey(rdmol)
            else:
                raise NotImplementedError('"{}" is not a valid output type'.format(otype))
        return odepic


    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of all cross references for a given chemical id (MNX) to other database id's
    #  Structure if the return: chebi_cid['88281']: 'MXM2323'
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers
    #TODO: save the self.deprecatedCID_cid to be used in case there rp_paths uses an old version of MNX
#    def _m_chebi_cid(self, cid_xref):
    @staticmethod
    def _m_chebi_cid(cid_xref):
        chebi_cid = {}
        for cid in cid_xref:
            if 'chebi' in cid_xref[cid]:
                for c in cid_xref[cid]['chebi']:
                    chebi_cid[c] = cid
        return chebi_cid

    ## Function to build the dictionnary to find the chemical id from inchikey
    #
    # @param cid_strc Dictionnary of chemical ID's to all the structure information associated with it
    # @return Dictionnary of InChIKey to chemical ID
    @staticmethod
    def _m_inchikey_cid(cid_strc):
        inchikey_cid = {}
        for cid in cid_strc:
            inchikey = cid_strc[cid]['inchikey']
            # This line is needed to put a value in 'inchikey', otherwise there are some problems in future strucutres
            if not inchikey: inchikey = 'NO_INCHIKEY'
            if not inchikey in inchikey_cid:
                inchikey_cid[inchikey] = []
            inchikey_cid[inchikey].append(cid)
        return inchikey_cid
