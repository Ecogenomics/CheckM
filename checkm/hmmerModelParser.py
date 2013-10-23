###############################################################################
#
# hmmer.py - parse a HMMER model file
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

class HmmModelError(Exception):
    pass

class HmmModel(object):
    """Parse a HMMER model file."""
    def __init__(self, keys, model=None):
        super(HmmModel, self).__init__()
        for key, value in keys.items():
            setattr(self, key, value)
        self.model = model

    def __len__(self):
        if self.leng is None:
            raise HmmModelError
        return self.leng

    def __str__(self):
        ret = str()

        print self.format
        ret += "NAME\t" + self.name + "\n"
        try:
            ret += "ACC\t" + self.acc + "\n"
        except AttributeError:
            pass
        try:
            ret += "DESC\t" + self.desc + "\n"
        except AttributeError:
            pass
        ret += "LENG\t"  + str(self.leng) + "\n"
        ret += "ALPH\t" + self.alph + "\n"
        try:
            if self.rf is True:
                ret += "RF\tyes\n"
            else:
                ret += "RF\tno\n"
        except AttributeError:
            pass
        try:
            if self.cs is True:
                ret += "CS\tyes\n"
            else:
                ret += "CS\tno\n"
        except AttributeError:
            pass
        try:
            if self.map is True:
                ret += "MAP\tyes\n"
            else:
                ret += "MAP\tno\n"
        except AttributeError:
            pass
        try:
            ret += "DATE\t" + self.date + "\n"
        except AttributeError:
            pass
        try:
            ret += "COM\t" + self.com + "\n"
        except AttributeError:
            pass
        try:
            ret += "NSEQ\t" + str(self.nseq) + "\n"
        except AttributeError:
            pass
        try:
            ret += "EFFN\t" + str(self.effn) + "\n"
        except AttributeError:
            pass
        try:
            ret += "CKSUM\t" + str(self.cksum) + "\n"
        except AttributeError:
            pass
        try:
            ret += "GA\t" + str(self.ga[0]) +" "+ str(self.ga[1]) + "\n"
        except AttributeError:
            pass
        try:
            ret += "TC\t" + str(self.tc[0]) +" "+ str(self.tc[1]) + "\n"
        except AttributeError:
            pass
        try:
            ret += "NC\t" + str(self.nc[0]) +" "+ str(self.nc[1]) + "\n"
        except AttributeError:
            pass
        try:
            ret += "STATS LOCAL MSV\t" + str(self.stats_local_msv[0]) +" "+ str(self.stats_local_msv[1]) + "\n"
        except AttributeError:
            pass
        try:
            ret += "STATS LOCAL VITERBI\t" + str(self.stats_local_viterbi[0]) +" "+ str(self.stats_local_viterbi[1]) + "\n"
        except AttributeError:
            pass
        try:
            ret += "STATS LOCAL FORWARD\t"+ str(self.stats_local_forward[0]) +" "+ str(self.stats_local_forward[1]) + "\n"
        except AttributeError:
            pass
        if self.model is not None:
            ret += self.model
        ret += "//\n"
        return ret


class HmmModelParser(object):
    """HmmModelParser holds a file object for a HMM model and a custom iterator
       for getting the values out"""
    def __init__(self, hmmfile):
        super(HmmModelParser, self).__init__()
        self.hmmfile = open(hmmfile)

    def parse(self):
        fields = []
        header_keys = dict()
        for current_line in self.hmmfile:
            # line should be: HMMER3/b [3.0b2 | June 2009]
            if current_line.startswith("HMMER"):
                header_keys['format'] = current_line.rstrip()

            elif current_line.startswith("HMM"):
                # begining of the model hmm
                # parsing not implemented at the moment - iterate through till
                # the end of this model
                model = ""
                for current_line in self.hmmfile:
                    if current_line.startswith("//"):
                        yield HmmModel(header_keys, model)
                        break
                    else:
                        model += current_line

            else:
                # header sections
                fields = current_line.rstrip().split(None, 1)
                if 2 != len(fields):
                    raise HmmModelError
                else:
                    # transform some data based on some of the header tags
                    if fields[0] == "LENG" or fields[0] == "NSEQ" or fields[0] == "CKSUM":
                        header_keys[fields[0].lower()] = int(fields[1])
                    elif fields[0] == "RF" or fields[0] == "CS" or fields[0] == "MAP":
                        if fields[1].lower() == "no":
                            header_keys[fields[0].lower()] = False
                        else:
                            header_keys[fields[0].lower()] = True
                    elif fields[0] == "EFFN":
                        header_keys[fields[0].lower()] = float(fields[1])
                    elif fields[0] == "GA" or fields[0] == "TC" or fields[0] == "NC":
                        params = fields[1].split()
                        if len(params) != 2:
                            raise HmmModelError
                        header_keys[fields[0].lower()] = (float(params[0].replace(';','')), float(params[1].replace(';','')))
                    elif fields[0] == "STATS":
                        params = fields[1].split()
                        if params[0] != "LOCAL":
                            raise HmmModelError
                        if params[1] == "MSV" or params[1] == "VITERBI" or params[1] == "FORWARD":
                            header_keys[(fields[0]+"_"+params[0]+"_"+params[1]).lower()] = (float(params[2]), float(params[3]))
                        else:
                            print("'"+params[1]+"'")
                            raise HmmModelError
                    else:
                        header_keys[fields[0].lower()] = fields[1]