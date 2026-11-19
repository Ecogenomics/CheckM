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
