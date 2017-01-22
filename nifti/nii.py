#!/usr/bin/env python
#
# Volconv - geometry-aware DICOM-to-NIfTI converter
# NIfTI module
#
# Copyright 2006-2017 Mark J White <mark@celos.net>
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
# See COPYING.txt and NOTICE.txt in the distribution for details.
#

from raw import *
import numpy

class NiiWriter(RawWriter):

    def __init__(self,filename,descrip=""):
        RawWriter.__init__(self,filename)

        self.offset = 352
        self.pixdim = [1.0, 1.0, 1.0] # not inc qfac
        self.qfac  = 1
        self.qform = 0
        self.sform = 0
        self.descrip = descrip
        self.qdata = [0.0,] * 6
        self.sdata = [0.0,] * 12
        self.one_padding = False

    def set_one_padding(self,value):
        self.one_padding = value

    def writeVal(self, addr, fmt, n, data, padding=None):

        # write the value
        self.fh.seek(addr)
        elts = len(data)
        str = self.order + fmt * elts
        bytes = struct.pack(str, *data)
        self.fh.write(bytes)

        if elts < n:

            # pad remaining length with zero values
            if padding is None:
                str = self.order + fmt * (n - elts)
                if fmt == 'c':
                    bytes = struct.pack(str, *('\x00' * (n - elts)))
                else:
                    bytes = struct.pack(str, *((0,) * (n - elts)))
                self.fh.write(bytes)

            # pad with given value repeatedly
            else:
                str = self.order + fmt * (n - elts)
                bytes = struct.pack(str, *((padding,) * (n - elts)))
                self.fh.write(bytes)


    def writeHeader(self):
        self.fh.seek(0)
        self.fh.write('\x00'*348)

        self.writeVal(0,   'i',  1, (348,)) # sizeof_hdr
        self.writeVal(40,  'h',  1, (numpy.rank(self.data),))

        if self.one_padding:
            self.writeVal(42,  'h',  7, numpy.shape(self.data), 1)
        else:
            self.writeVal(42,  'h',  7, numpy.shape(self.data))

        self.writeVal(70,  'h',  1, (self.type,))
        self.writeVal(72,  'h',  1, (NiftiType.Map[self.type][3],))

        if self.one_padding:
            self.writeVal(80,  'f',  7, self.pixdim, 1.0)
        else:
            self.writeVal(80,  'f',  7, self.pixdim)

        self.writeVal(108, 'f',  7, (352.0,)) # vox_offset

        if self.one_padding:
            self.writeVal(112, 'f',  1, (1.0,)) # scl_slope
            self.writeVal(116, 'f',  1, (0.0,)) # scl_inter

        self.writeVal(123, 'b',  1, (10,)) # xyzt_units: 10 = mm + sec

        self.writeVal(344, 'c',  4, "n+1\x00") # magic
        self.writeVal(148, 'c', 80, self.descrip)

        self.writeVal(252, 'h',  1, (self.qform,))
        self.writeVal(254, 'h',  1, (self.sform,))

        self.writeVal(256, 'f',  6, self.qdata)
        self.writeVal(76,  'f',  1, (self.qfac,))

    def write(self):
        self.writeHeader()
        RawWriter.write(self)
