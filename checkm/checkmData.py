###############################################################################
#
# checkmData.py - database management utilities
#
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import sys
import logging

try: # setuptools<81
    from pkg_resources import resource_filename
    def get_resource(resource, pkg='checkm'):
        return os.path.abspath(resource_filename(pkg, resource))

except(ModuleNotFoundError): # setuptools>=81
    from importlib import resources
    def get_resource(resource, pkg='checkm'):
        return resources.files(pkg) / resource

import json
from collections import namedtuple

import checkm.manifestManager as mm


class DBConfig(object):
    """CheckM uses setuptools to distribute a file called "DATA_CONFIG" which is placed with
    the other files during installation. This file stores information about where data is
    stored locally and where to look for remote updates. This class essentially exposes
    the DATA_CONFIG file as an object.
    """

    def __init__(self, configFile):
        self.logger = logging.getLogger('timestamp')
        self.configFile = configFile
        self.values = self.getConfig()

# -----------------------------------------------------------------------------
# Read and write local config file

    def getConfig(self):
        """Get a listing of the versions of files in the local data config"""
        try:
            with open(self.configFile, 'r') as local_config:
                # config is a one line file
                for line in local_config:
                    return json.loads(line)

        except Exception:
            print("There seems to be a problem with loading the CheckM config file")
            print("Please check the permissions / existence / contents of:")
            print(self.configFile)
            raise

        return {}

    def setConfig(self):
        """Update the local config to reflect any changes made"""
        if self.checkPermissions():
            with open(self.configFile, 'w') as config_fh:
                config_fh.write(json.dumps(self.values))

# -----------------------------------------------------------------------------
# Housekeeping

    def checkPermissions(self):
        """Work out if we have permission to write to the CheckM config before attempting to make changes"""
        try:
            with open(self.configFile, 'a'):
                pass
        except IOError as e:
            print("You do not seem to have permission to edit the CheckM config file")
            print("located at %s" % self.configFile)
            print("Please try again with updated privileges. Error was:\n")
            print(e)
            return False
        return True


class DBManager(mm.ManifestManager):

    """Manage all aspects of data location and version control."""

    def __init__(self, set_path=None, configFile=None):
        mm.ManifestManager.__init__(self)

        self.logger = logging.getLogger('timestamp')

        # check if CheckM data path is indicated by environment variable
        if 'CHECKM_DATA_PATH' in os.environ:
            data_root = os.environ['CHECKM_DATA_PATH']

            Config = namedtuple('Config', 'values')
            self.config = Config(values={'dataRoot': data_root})

            if set_path:
                self.logger.warning(
                    'The CHECKM_DATA_PATH environment variable is set (%s), so the data root in the CheckM '
                    'DATA_CONFIG file is not used and has not been changed. Unset CHECKM_DATA_PATH to use '
                    "'checkm data setRoot', or set CHECKM_DATA_PATH to the desired directory." % data_root)
        else:
            # read data path from DATA_CONFIG file
            if configFile is None:
                configFile = get_resource('DATA_CONFIG')

            self.config = DBConfig(configFile=configFile)
            self.type = self.config.values["manifestType"]

            if set_path:
                if not self.setRoot(set_path):
                    print("Unable to set the CheckM data folder to: %s" % set_path)
                    sys.exit(1)

            # check that the data root is legit
            manifestFile = os.path.join(
                self.config.values["dataRoot"], mm.__MANIFEST__)
            if not os.path.exists(self.config.values["dataRoot"]) or not os.path.exists(manifestFile):
                self.config.values["dataRoot"] = ""

            if self.config.values["dataRoot"] == "":
                # no data folder set.
                print("It seems that the CheckM data folder has not been set yet or has been removed. Please run 'checkm data setRoot'.")
                if not self.setRoot():
                    print("CheckM cannot run without a valid data folder.")

    def setRoot(self, path=None):
        """Set the data folder"""
        # check to see we have permission to update the data root
        if not self.config.checkPermissions():
            return None

        # validate the supplied path or prompt for a new one
        path = self.confirmPath(path=path)
        if path is None:
            # The user is not interested
            return None

        # path should be set, exist and be writable
        self.config.values["dataRoot"] = path

        # save the new path
        self.config.setConfig()

        return path

    def confirmPath(self, path=None):
        """Check that the data path exists and is writable, creating it if necessary.

        Returns the absolute path, or None if the path could not be created.
        """
        if not path:
            path = os.path.join(os.path.expanduser("~"), ".checkm")

        if path.upper() == "ABORT":
            # user has given up
            return None

        path = os.path.abspath(os.path.expanduser(path))

        print("")
        if os.path.exists(path):
            # path exists
            if os.access(path, os.W_OK):
                # path is writable
                print(
                    "Path [%s] exists and you have permission to write to this folder." % path)
            else:
                print(
                    "Path [%s] exists but you do not have permission to write to this folder." % path)
                sys.exit(1)
        else:
            # path does not exist, try to make it
            print(
                "Path [%s] does not exist so I will attempt to create it" % path)
            try:
                self.makeSurePathExists(path)
                print(
                    "Path [%s] has been created and you have permission to write to this folder." % path)
            except Exception as e:
                print("Unable to make the folder; error was: %s" % e)
                return None

        # (re)make the manifest file
        print("(re) creating manifest file (please be patient).")
        self.createManifest(path, self.config.values["localManifestName"])

        return path
