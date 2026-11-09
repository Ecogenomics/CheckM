###############################################################################
#
# progress.py - progress bars for long-running steps
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

import sys
import logging

from tqdm import tqdm


def progressEnabled(logger=None):
    """Check if progress should be reported to the console.

    Progress is reported when the 'timestamp' logger reports INFO messages to the
    console, i.e., unless the logger's level is above INFO or the console output
    has been suppressed (e.g., with --quiet).
    """
    logger = logger or logging.getLogger('timestamp')
    if logger.getEffectiveLevel() > logging.INFO:
        return False

    for handler in logger.handlers:
        if (isinstance(handler, logging.StreamHandler)
                and not isinstance(handler, logging.FileHandler)
                and handler.level > logging.INFO):
            return False

    return True


def progressBar(total, desc, unit, enabled=True):
    """Create a progress bar written to stderr (e.g., '    Processing bins: 100%|...| 8/8 [...]').

    Parameters
    ----------
    total : int
        Number of items to process.
    desc : str
        Description of the step (e.g., 'Processing bins').
    unit : str
        Name of a single item (e.g., 'bin').
    enabled : bool
        Set to False to never show the progress bar.

    The progress bar is disabled if progress should not be reported to the
    console (see progressEnabled), but can still be updated and closed.
    """
    return tqdm(total=total,
                desc='    ' + desc,
                unit=unit,
                file=sys.stderr,
                leave=True,
                disable=not (enabled and progressEnabled()))
