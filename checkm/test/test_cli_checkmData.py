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

"""Tests for setting the CheckM data root (checkm data setRoot).

All tests work on a temporary copy of DATA_CONFIG and a temporary HOME, so
neither the packaged DATA_CONFIG file nor ~/.checkm is modified.
"""

import os
import json
import shutil
import tempfile
import unittest
from unittest import mock

from checkm.checkmData import DBManager, DBConfig
from checkm.manifestManager import ManifestManager


class TestDBManager(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()
        self.home = os.path.join(self.tmpDir, 'home')
        os.makedirs(self.home)

        self.configFile = os.path.join(self.tmpDir, 'DATA_CONFIG')
        with open(self.configFile, 'w') as f:
            f.write(json.dumps({"dataRoot": "",
                                "remoteManifestURL": "https://example.org/",
                                "manifestType": "CheckM",
                                "localManifestName": ".dmanifest",
                                "remoteManifestName": ".dmanifest"}))

        self.dataDir = os.path.join(self.tmpDir, 'data')
        os.makedirs(os.path.join(self.dataDir, 'hmms'))
        with open(os.path.join(self.dataDir, 'hmms', 'a.hmm'), 'w') as f:
            f.write('HMM')
        with open(os.path.join(self.dataDir, 'b.tsv'), 'w') as f:
            f.write('b')

        # never touch the real environment, HOME or DATA_CONFIG
        env = {k: v for k, v in os.environ.items() if k != 'CHECKM_DATA_PATH'}
        env['HOME'] = self.home
        self.envPatch = mock.patch.dict(os.environ, env, clear=True)
        self.envPatch.start()

    def tearDown(self):
        self.envPatch.stop()
        shutil.rmtree(self.tmpDir)

    def readConfig(self):
        with open(self.configFile) as f:
            return json.loads(f.readline())

    def testSetRoot(self):
        dbm = DBManager(set_path=self.dataDir, configFile=self.configFile)
        self.assertEqual(dbm.config.values['dataRoot'], self.dataDir)
        self.assertEqual(self.readConfig()['dataRoot'], self.dataDir)

        # every file is listed exactly once in the manifest
        with open(os.path.join(self.dataDir, '.dmanifest')) as f:
            entries = [line.split('\t')[0] for line in f if not line.startswith('#')]
        self.assertEqual(len(entries), len(set(entries)))
        self.assertEqual(sorted(entries), sorted(['b.tsv', 'hmms', os.path.join('hmms', 'a.hmm')]))

        self.assertFalse(os.path.exists(os.path.join(self.home, '.checkm')))

    def testSetRootUncreatablePath(self):
        notADir = os.path.join(self.tmpDir, 'file')
        with open(notADir, 'w') as f:
            f.write('x')

        with self.assertRaises(SystemExit) as cm:
            DBManager(set_path=os.path.join(notADir, 'data'), configFile=self.configFile)
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(self.readConfig()['dataRoot'], '')
        self.assertFalse(os.path.exists(os.path.join(self.home, '.checkm')))

    def testConfirmPathUncreatablePath(self):
        DBManager(set_path=self.dataDir, configFile=self.configFile)

        notADir = os.path.join(self.tmpDir, 'file')
        with open(notADir, 'w') as f:
            f.write('x')

        dbm = DBManager(configFile=self.configFile)
        self.assertIsNone(dbm.confirmPath(os.path.join(notADir, 'data')))

    def testEnvironmentVariableOverridesConfig(self):
        os.environ['CHECKM_DATA_PATH'] = self.dataDir
        with self.assertLogs('timestamp', level='WARNING') as logs:
            dbm = DBManager(set_path=os.path.join(self.tmpDir, 'other'), configFile=self.configFile)
        self.assertIn('CHECKM_DATA_PATH', '\n'.join(logs.output))
        self.assertEqual(dbm.config.values['dataRoot'], self.dataDir)
        self.assertEqual(self.readConfig()['dataRoot'], '')
        self.assertFalse(os.path.exists(os.path.join(self.tmpDir, 'other')))

    def testCheckPermissions(self):
        config = DBConfig(self.configFile)
        self.assertTrue(config.checkPermissions())

        config.configFile = os.path.join(self.tmpDir, 'missing_dir', 'DATA_CONFIG')
        self.assertFalse(config.checkPermissions())


class TestManifestManager(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.tmpDir, 'sub'))
        with open(os.path.join(self.tmpDir, 'sub', 'x.txt'), 'w') as f:
            f.write('x')

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def testRepeatedCreateManifest(self):
        mm = ManifestManager()
        mm.createManifest(self.tmpDir)
        mm.createManifest(self.tmpDir)

        with open(os.path.join(self.tmpDir, '.dmanifest')) as f:
            entries = [line.split('\t')[0] for line in f if not line.startswith('#')]
        self.assertEqual(sorted(entries), sorted(['sub', os.path.join('sub', 'x.txt')]))


if __name__ == "__main__":
    unittest.main()
