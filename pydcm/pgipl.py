#!/usr/bin/env python
#
# Volconv - geometry-aware DICOM-to-NIfTI converter
# Slow, simple GIPL writer
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


import struct
import sys
import numpy

class GiplType:

    # Typecodes
    Binary    = 1
    Char      = 7
    UChar     = 8
    Short     = 15
    UShort    = 16
    Int       = 31
    UInt      = 32
    Float     = 64
    Double    = 65
    CShort    = 144
    CInt      = 160
    CFloat    = 192
    CDouble   = 193

    # Magic numbers
    Magic     = 719555000
    Magic2    = 4026526128
    MagicEnd  = 3096044330
    Magic2End = 2968125423
    MagicExt  = 815

class GiplFile:

    def __init__(self, filename, mode):
        self.fn = filename
        self.mode = mode
        if self.mode == "r":
            self.fh = file(filename,"rb")
        elif self.mode == "w":
            self.fh = file(filename,"wb")
        else:
            raise RuntimeError("Unknown mode")

        self.IS = (1,1,1,1)
        self.DT = GiplType.Char
        self.VS = (1., 1., 1., 1.)
        self.DE = ""
        self.MA = (1., 0., 0., 0.,
                        0., 1., 0., 0.,
                        0., 0., 1., 0.)
        self.RL = 0.
        self.RH = 0.
        self.OF = (0., 0., 0., 0.)

        # if sys.byteorder == 'big':
        #     self.order = ">"
        # else:
        #     self.order = "<"
        self.order = ">"

        self.ext = True

    def readVal(self, addr, fmt, n):
        self.fh.seek(addr)
        str = self.order + fmt * n
        bytes = struct.calcsize(str)
        string = self.fh.read(bytes)
        return struct.unpack(str, string)

    def writeVal(self, addr, fmt, n, data, seq):
        self.fh.seek(addr)
        str = self.order + fmt * n
        if seq:
            data = struct.pack(str, *data)
        else:
            data = struct.pack(str, data)
        self.fh.write(data)
    
    def zeroHeader(self):
        self.fh.seek(0)
        self.fh.write("\x00"*256)

    def dump(self):
        print "File:                ", self.fn
        print "  Byte order:        ", self.order
        print "  IS (image size):   ", self.IS
        print "  DR (data type):    ", self.DT
        print "  VS (voxel size):   ", self.VS
        print "  DE (description):  ", self.DE

        if (self.ext):
            print "  MA (matrix):       ", self.MA
            print "  OF (offset origin):", self.OF
            print "  RL (minimum):      ", self.RL
            print "  RH (maximum):      ", self.RH

    def naType(self):
        if self.DT == GiplType.Binary:
            return "Bool"
        elif self.DT == GiplType.Char:
            return "Int8"
        elif self.DT == GiplType.UChar:
            return "UInt8"
        elif self.DT == GiplType.Short:
            return "Int16"
        elif self.DT == GiplType.UShort:
            return "UInt16"
        elif self.DT == GiplType.Int:
            return "UInt32"
        elif self.DT == GiplType.UInt:
            return "UInt64"
        elif self.DT == GiplType.Float:
            return "Float32"
        elif self.DT == GiplType.Double:
            return "Float64"
        elif self.DT == GiplType.CFloat:
            return "Complex32"
        elif self.DT == GiplType.CDouble:
            return "Complex64"
        elif (self.DT == GiplType.CInt
                or self.DT == GiplType.CShort):
            raise RuntimeError("No complex int or short with numarray yet")
        else:
            raise RuntimeError("Unknown type")

    def packType(self):
        if self.DT == GiplType.Binary:
            return "c"
        elif self.DT == GiplType.Char:
            return "b"
        elif self.DT == GiplType.UChar:
            return "B"
        elif self.DT == GiplType.Short:
            return "h"
        elif self.DT == GiplType.UShort:
            return "H"
        elif self.DT == GiplType.Int:
            return "i"
        elif self.DT == GiplType.UInt:
            return "H"
        elif self.DT == GiplType.Float:
            return "f"
        elif self.DT == GiplType.Double:
            return "d"
        elif self.DT == GiplType.CFloat:
            return "ff"
        elif self.DT == GiplType.CDouble:
            return "dd"
        elif (self.DT == GiplType.CInt
                or self.DT == GiplType.CShort):
            raise RuntimeError("No complex int or short with numarray yet")
        else:
            raise RuntimeError("Unknown type")

    def isComplex(self):
        if (self.DT == GiplType.CFloat 
                or self.DT == GiplType.CDouble
                or self.DT == GiplType.CInt
                or self.DT == GiplType.CShort):
            return True
        else:
            return False

    def alloc(self):
        self.data = numarray.zeros(self.IS, self.naType())

    def readData(self):
        self.alloc()
        linear = numarray.reshape(self.data, 
                (self.IS[0]*self.IS[1]*self.IS[2], self.IS[3]))
        self.fh.seek(256)

        planen = self.IS[0] * self.IS[1] * self.IS[2]
        str = self.order + self.packType() * planen
        bytes = struct.calcsize(str)

        for i in range(0, self.IS[3]):
            tmp = self.fh.read(bytes)
            d = struct.unpack(str, tmp)
            # XXX deal with complex
            linear[:,i] = d

    def writeData(self):
        if self.mode == "r":
            raise RuntimeError("Called writeData() with read mode")
        
        linear = numarray.reshape(self.data, 
                (self.IS[0]*self.IS[1]*self.IS[2], self.IS[3]))
        self.fh.seek(256)

        planen = self.IS[0] * self.IS[1] * self.IS[2]
        str = self.order + self.packType() * planen
        bytes = struct.calcsize(str)

        for i in range(0, self.IS[3]):
            bytes = struct.pack(str, *linear[:,i])
            self.fh.write(bytes)

    def read(self):
        self.validFormat()
        
        self.IS = self.readVal(0,  "H", 4)
        self.DT = self.readVal(8,  "H", 1)[0]
        self.VS = self.readVal(10, "f", 4)
        self.DE = self.readVal(26, "c", 80)[0].rstrip("\x00")

        if self.readVal(244, "H", 1)[0] == GiplType.MagicExt:
            self.RL = self.readVal(188,"d", 1)[0]
            self.RH = self.readVal(196,"d", 1)[0]
            self.MA = self.readVal(106,"f", 12)
            self.OF = self.readVal(204,"d", 4)
            self.ext = True
        else:
            self.ext = False

        self.readData()

    def write(self):
        if self.mode == "r":
            raise RuntimeError("Called write() with read mode")

        # clear header
        self.zeroHeader()

        # write data block
        self.writeVal(0,   "H", 4, self.IS, 1)
        self.writeVal(8,   "H", 1, self.DT, 0)
        self.writeVal(10,  "f", 4, self.VS, 1)
        self.writeVal(26,  "c",80, self.DE + "\x00" * (80-len(self.DE)), 1)
        self.writeVal(106, "f",12, self.MA, 1)
        self.writeVal(188, "d", 1, self.RL, 0)
        self.writeVal(196, "d", 1, self.RH, 0)
        self.writeVal(204, "d", 4, self.OF, 1)

        # write magic flags
        self.writeVal(252, "I", 1, GiplType.Magic, 0)
        if self.ext:
            self.writeVal(244, "I", 1, GiplType.MagicExt, 0)

        self.writeData()

    def validFormat(self):
        self.fh.seek(252)
        magic = struct.unpack(self.order + "I", self.fh.read(4))[0]
        if magic == GiplType.Magic or magic == GiplType.Magic2:
            return True
        if magic == GiplType.MagicEnd or magic == GiplType.Magic2End:
            if sys.byteorder == 'big':
                self.order = "<"
            else:
                self.order = ">"
            return True
        else:
            raise RuntimeError("Not a GIPL file")
