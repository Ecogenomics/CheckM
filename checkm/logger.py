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
import ntpath


def logger_setup(log_dir, log_file, program_name, version, silent):
    """Set up the 'timestamp' logger.

    The 'timestamp' logger prefixes each message with a timestamp and
    writes to stdout and, when log_dir is not None, to a log file in
    log_dir. If the log file cannot be created, a warning is reported
    and messages are only written to stdout. Calling this function
    again replaces any handlers added by a previous call, so messages
    are never reported twice.

    Parameters
    ----------
    log_dir : str
        Output directory for log file.
    log_file : str
        Desired name of log file.
    program_name : str
        Name of program.
    version : str
        Program version number.
    silent : boolean
        Flag indicating if output to stdout should be suppressed.
    """

    # set up general properties of logger
    timestamp_logger = logging.getLogger('timestamp')
    timestamp_logger.setLevel(logging.DEBUG)
    log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                   datefmt="%Y-%m-%d %H:%M:%S")

    # remove handlers from any previous call
    for handler in list(timestamp_logger.handlers):
        timestamp_logger.removeHandler(handler)
        handler.close()

    # set up logging to console
    timestamp_stream_logger = logging.StreamHandler(sys.stdout)
    timestamp_stream_logger.setFormatter(log_format)
    timestamp_logger.addHandler(timestamp_stream_logger)

    if silent:
        timestamp_stream_logger.setLevel(logging.ERROR)

    if log_dir:
        try:
            if not os.path.exists(log_dir):
                os.makedirs(log_dir)
            timestamp_file_logger = logging.FileHandler(os.path.join(log_dir, log_file), 'a')
            timestamp_file_logger.setFormatter(log_format)
            timestamp_logger.addHandler(timestamp_file_logger)
        except OSError as e:
            timestamp_logger.warning('Unable to write log file to %s: %s' % (log_dir, e))

    timestamp_logger.info('%s v%s' % (program_name, version))
    timestamp_logger.info(ntpath.basename(sys.argv[0]) + ' ' + ' '.join(sys.argv[1:]))
