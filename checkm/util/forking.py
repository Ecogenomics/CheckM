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

import gc
import sys
import multiprocessing as mp
from contextlib import contextmanager


@contextmanager
def frozenHeap():
    """Start worker processes without them later copying the memory of this process.

    Forked workers share the memory of their parent until they write to it. A
    full garbage collection in a worker writes to the header of every object it
    inherited, which leaves each worker with a private copy of the parent's
    data. Objects frozen with gc.freeze() are ignored by the collector, so
    worker processes started within this context keep sharing them.
    """

    gc.freeze()
    try:
        yield
    finally:
        gc.unfreeze()


def runInChild(func, *args, **kwargs):
    """Run a function in a child process that exits once the function returns.

    HMMER is run within the Python process (through pyhmmer). The memory of a
    worker process that runs HMMER many times grows well beyond that needed
    for a single run, and the process writes to memory it shares with its
    parent, which is then copied. Running each task in its own short-lived
    process returns this memory to the system, as when the HMMER programs
    were run as separate processes. Errors are reported by the function,
    which must exit with a non-zero status on failure; this process then
    exits with status 1. The function is run directly if processes can't be
    forked.
    """
    if 'fork' not in mp.get_all_start_methods():
        func(*args, **kwargs)
        return

    proc = mp.get_context('fork').Process(target=func, args=args, kwargs=kwargs)
    proc.start()
    proc.join()
    if proc.exitcode != 0:
        sys.exit(1)
