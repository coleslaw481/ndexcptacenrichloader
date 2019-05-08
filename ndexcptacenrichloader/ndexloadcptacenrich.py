#! /usr/bin/env python

import argparse
import sys
import os
import logging
import json
from logging import config
from ndexutil.config import NDExUtilConfig
import ndex2
from ndex2.client import Ndex2
import ndexcptacenrichloader

logger = logging.getLogger(__name__)

TSV2NICECXMODULE = 'ndexutil.tsv.tsv2nicecx2'

LOG_FORMAT = "%(asctime)-15s %(levelname)s %(relativeCreated)dms " \
             "%(filename)s::%(funcName)s():%(lineno)d %(message)s"


TYPE_MAP = {'GeneProduct': 'protein',
            'Protein': 'protein',
            'Group': 'complex'}

STYLE = 'style.cx'
"""
Name of file containing CX with style
stored within this package
"""


def get_package_dir():
    """
    Gets directory where package is installed
    :return:
    """
    return os.path.dirname(ndexcptacenrichloader.__file__)


def get_style():
    """
    Gets the style stored with this package
    :return: path to file
    :rtype: string
    """
    return os.path.join(get_package_dir(), STYLE)


def _parse_arguments(desc, args):
    """
    Parses command line arguments
    :param desc:
    :param args:
    :return:
    """
    help_fm = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_fm)
    parser.add_argument('--profile', help='Profile in configuration '
                                          'file to use to load '
                                          'NDEx credentials which means'
                                          'configuration under [XXX] will be'
                                          'used '
                                          '(default '
                                          'ndexcptacenrichloader)',
                        default='ndexcptacenrichloader')
    parser.add_argument('--logconf', default=None,
                        help='Path to python logging configuration file in '
                             'this format: https://docs.python.org/3/library/'
                             'logging.config.html#logging-config-fileformat '
                             'Setting this overrides -v parameter which uses '
                             ' default logger. (default None)')

    parser.add_argument('--conf', help='Configuration file to load '
                                       '(default ~/' +
                                       NDExUtilConfig.CONFIG_FILE)
    parser.add_argument('--style',
                        help='Path to NDEx CX file to use for styling '
                             'networks', default=get_style())
    parser.add_argument('--verbose', '-v', action='count', default=0,
                        help='Increases verbosity of logger to standard '
                             'error for log messages in this module and'
                             'in ' + TSV2NICECXMODULE + '. Messages are '
                             'output at these python logging levels '
                             '-v = ERROR, -vv = WARNING, -vvv = INFO, '
                             '-vvvv = DEBUG, -vvvvv = NOTSET (default no '
                             'logging)')
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' +
                                 ndexcptacenrichloader.__version__))

    return parser.parse_args(args)


def _setup_logging(args):
    """
    Sets up logging based on parsed command line arguments.
    If args.logconf is set use that configuration otherwise look
    at args.verbose and set logging for this module and the one
    in ndexutil specified by TSV2NICECXMODULE constant
    :param args: parsed command line arguments from argparse
    :raises AttributeError: If args is None or args.logconf is None
    :return: None
    """

    if args.logconf is None:
        level = (50 - (10 * args.verbose))
        logging.basicConfig(format=LOG_FORMAT,
                            level=level)
        logging.getLogger(TSV2NICECXMODULE).setLevel(level)
        logger.setLevel(level)
        return

    # logconf was set use that file
    logging.config.fileConfig(args.logconf,
                              disable_existing_loggers=False)


class NDExCPTACLoader(object):
    """
    Class to load content
    """

    SOURCE_PREFIX='source_'
    NETWORKSET='networkset'

    def __init__(self, args):
        """

        :param args:
        """
        self._args = args
        self._conf_file = args.conf
        self._profile = args.profile
        self._ndex = None
        self._sourcendex = None
        self._user = None
        self._pass = None
        self._server = None
        self._sourceuser = None
        self._sourcepass = None
        self._sourceserver = None
        self._sourcenetworkset = None
        self._template = None

    def _parse_config(self):
            """
            Parses config
            :return:
            """
            ncon = NDExUtilConfig(conf_file=self._conf_file)
            con = ncon.get_config()
            self._user = con.get(self._profile, NDExUtilConfig.USER)
            self._pass = con.get(self._profile, NDExUtilConfig.PASSWORD)
            self._server = con.get(self._profile, NDExUtilConfig.SERVER)
            self._sourceuser = con.get(self._profile,
                                       NDExCPTACLoader.SOURCE_PREFIX +
                                       NDExUtilConfig.USER)
            self._sourcepass = con.get(self._profile,
                                       NDExCPTACLoader.SOURCE_PREFIX +
                                       NDExUtilConfig.PASSWORD)
            self._sourceserver = con.get(self._profile,
                                         NDExCPTACLoader.SOURCE_PREFIX +
                                         NDExUtilConfig.SERVER)
            self._sourcenetworkset = con.get(self._profile,
                                         NDExCPTACLoader.SOURCE_PREFIX +
                                         NDExCPTACLoader.NETWORKSET)

    def _get_user_agent(self):
        """

        :return:
        """
        return 'cptac/' + self._args.version

    def _create_ndex_connection(self):
        """
        creates connection to ndex
        :return:
        """
        if self._ndex is None:
            self._ndex = Ndex2(host=self._server, username=self._user,
                               password=self._pass,
                               user_agent=self._get_user_agent())

    def _create_sourcendex_connection(self):
        """
        creates connection to ndex
        :return:
        """
        if self._sourcendex is None:
            self._sourcendex = Ndex2(host=self._sourceserver, username=self._sourceuser,
                                     password=self._sourcepass,
                                     user_agent=self._get_user_agent())

    def _load_style_template(self):
        """
        Loads the CX network specified by self._args.style into self._template
        :return:
        """
        self._template = ndex2.create_nice_cx_from_file(os.path.abspath(self._args.style))

    def _load_network_summaries_for_user(self):
        """
        Gets a dictionary of all networks for user account
        <network name upper cased> => <NDEx UUID>
        :return: dict
        """
        net_summaries = self._ndex.get_network_summaries_for_user(self._user)
        self._net_summaries = {}
        for nk in net_summaries:
            if nk.get('name') is not None:
                self._net_summaries[nk.get('name').upper()] = nk.get('externalId')

    def _get_list_of_networkids(self):
        """
        Queries source ndex client with networkset passed in configuration to get
        list of NDEx network UUIDs
        :return:
        """
        res = self._sourcendex.get_network_set(self._sourcenetworkset)
        if res is None:
            logger.error('No networks found')
            return None
        if 'networks' not in res:
            logger.error('No networks entry in result: ' + str(res))

        return res['networks']

    def _remap_raw_type_new_normalized_type(self, raw_type):
        """
        Uses map at top of this module to map raw_type
        to new type. If no match, raw_type is passed on
        :param raw_type:
        :return:
        """
        if raw_type is None:
            return None
        if raw_type in TYPE_MAP:
            return TYPE_MAP[raw_type]
        return raw_type

    def _process_network_by_id(self, networkid):
        """
        Processes network by id
        :param networkid:
        :return:
        """
        if networkid is None:
            logger.error('network id is None')
            return
        network = ndex2.create_nice_cx_from_server(self._sourceserver, username=self._sourceuser,
                                                   password=self._sourcepass,
                                                   uuid=networkid)

        logger.info('NETWORK: ' + network.get_name())
        for id, node in network.get_nodes():

            raw_type = network.get_node_attribute(id, 'WP.type')
            if raw_type is None:
                continue
            raw_type = raw_type['v']

            # some nodes have empty string for name which currently screws up enrichment
            # so going to just set the type to something else for these right now
            if 'n' not in node or node['n'] is None or len(node['n']) == 0:
                node['n'] = 'unset'
                network.add_node_attribute(property_of=id, name='type', values='unsetname' + raw_type, overwrite=True)
            else:
                network.add_node_attribute(property_of=id, name='type',
                                           values=self._remap_raw_type_new_normalized_type(raw_type),
                                           overwrite=True)

        # apply style to network
        network.apply_style_from_network(self._template)

        network_update_key = self._net_summaries.get(network.get_name().upper())

        if network_update_key is not None:
            return network.update_to(network_update_key, self._server, self._user, self._pass,
                                     user_agent=self._get_user_agent())
        else:
            upload_message = network.upload_to(self._server, self._user,
                                               self._pass,
                                               user_agent=self._get_user_agent())
        return upload_message

    def run(self):
        """
        Runs content loading for NDEx Cancer Hallmark networks from CPTAC at WikiPathways Enrichment Loader
        :param theargs:
        :return:
        """
        self._parse_config()
        self._create_ndex_connection()
        self._create_sourcendex_connection()
        logger.debug('Parsed config: ' + self._user)
        self._load_network_summaries_for_user()
        self._load_style_template()

        networklist = self._get_list_of_networkids()
        if networklist is None:
            logger.error('No networks')
            return 1

        for networkid in networklist:
            logger.debug('Processing network: ' + networkid)
            self._process_network_by_id(networkid)

        return 0


def main(args):
    """
    Main entry point for program
    :param args:
    :return:
    """
    desc = """
    Version {version}

    Loads NDEx Cancer Hallmark networks from CPTAC at WikiPathways Enrichment Loader data into NDEx (http://ndexbio.org).
    
    To connect to NDEx server a configuration file must be passed
    into --conf parameter. If --conf is unset the configuration 
    the path ~/{confname} is examined. 
         
    The configuration file should be formatted as follows:
         
    [<value in --profile (default ncipid)>]
         
    {user} = <NDEx username>
    {password} = <NDEx password>
    {server} = <NDEx server(omit http) ie public.ndexbio.org>
    {sourceuser} = <NDEx username>
    {sourcepass} = <NDEx password>
    {sourceserver} = <NDEx server (omit http) ie public.ndexbio.org>
    {sourcenetworksetid} = <networkset id containing networks to import>
    
    
    """.format(confname=NDExUtilConfig.CONFIG_FILE,
               user=NDExUtilConfig.USER,
               password=NDExUtilConfig.PASSWORD,
               server=NDExUtilConfig.SERVER,
               sourceuser=NDExCPTACLoader.SOURCE_PREFIX + NDExUtilConfig.USER,
               sourcepass=NDExCPTACLoader.SOURCE_PREFIX + NDExUtilConfig.PASSWORD,
               sourceserver=NDExCPTACLoader.SOURCE_PREFIX + NDExUtilConfig.SERVER,
               sourcenetworksetid=NDExCPTACLoader.SOURCE_PREFIX +NDExCPTACLoader.NETWORKSET,
               version=ndexcptacenrichloader.__version__)
    theargs = _parse_arguments(desc, args[1:])
    theargs.program = args[0]
    theargs.version = ndexcptacenrichloader.__version__

    try:
        _setup_logging(theargs)
        loader = NDExCPTACLoader(theargs)
        return loader.run()
    except Exception as e:
        logger.exception('Caught exception')
        return 2
    finally:
        logging.shutdown()


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
