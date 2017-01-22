#!/usr/bin/env python
#
# Volconv - geometry-aware DICOM-to-NIfTI converter
# GIPL header handling module
#
# Copyright 2006-2016 Mark J White <mark@celos.net>
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

class GiplType:

    # Typecodes
    Map = {
        NiftiType.Bool:     1,
        NiftiType.Int8:     7,
        NiftiType.UInt8:    8,
        NiftiType.Int16:    15,
        NiftiType.UInt16:   16,
        NiftiType.Int32:    31,
        NiftiType.UInt32:   32,
        NiftiType.Float32:  64,
        NiftiType.Float64:  65,
        NiftiType.Complex64:192,
        NiftiType.Complex128:193
    }

    # Magic numbers
    Magic     = 719555000
    Magic2    = 4026526128
    MagicEnd  = 3096044330
    Magic2End = 2968125423
    MagicExt  = 815


class GiplWriter(RawWriter):

    def __init__(self,filename):
        RawWriter.__init__(self,filename)

        self.offset = 256
        self.pixdim = [1.0, 1.0, 1.0] # not inc qfac
        self.qfac  = 1
        self.qform = 0
        self.sform = 0
        self.descrip = ""
        self.qdata = [0.0,] * 6
        self.sdata = [0.0,] * 12
        self.order = '>'

    def writeVal(self, addr, fmt, n, data):

        # write the value
        self.fh.seek(addr)
        elts = len(data)
        str = self.order + fmt * elts
        bytes = struct.pack(str, *data)
        self.fh.write(bytes)

        # pad remaining length with zero values
        if elts < n:
            str = self.order + fmt * (n - elts)
            if fmt == 'c':
                bytes = struct.pack(str, *('\x00' * (n - elts)))
            else:
                bytes = struct.pack(str, *((0,) * (n - elts)))
            self.fh.write(bytes)

    def writeHeader(self):
        self.fh.seek(0)
        self.fh.write('\x00'*256)

        my_pixdim = [1.0,1.0,1.0,1.0]
        my_pixdim[0:len(self.pixdim)] = self.pixdim[0:4]

        # write data block
        self.writeVal(0,   "H", 4, numpy.shape(self.data)) # dim
        self.writeVal(8,   "H", 1, (GiplType.Map[self.type],)) # type
        self.writeVal(10,  "f", 4, my_pixdim) # pixdim
        self.writeVal(26,  "c",80, self.descrip) # description
        self.writeVal(106, "f",12, (0,)) # matrix (ignore for now)
        self.writeVal(188, "d", 1, (0,)) # low
        self.writeVal(196, "d", 1, (0,)) # high
        self.writeVal(204, "d", 4, (0,)) # origin

        # write magic flags
        self.writeVal(252, "I", 1, (GiplType.Magic,))
        # if self.ext:
        #     self.writeVal(244, "I", 1, (GiplType.MagicExt,))

    def write(self):
        self.writeHeader()
        RawWriter.write(self)
