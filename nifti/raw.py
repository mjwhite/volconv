#!/usr/bin/env python
#
# Volconv - geometry-aware DICOM-to-NIfTI converter
# Raw writer (primarily for NIfTI data, hence using NIfTI type system)
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
import gzip

class NiftiType:

    # complete set of allowed Nifti-1.1 types

    Bool       =    1 
    UInt8      =    2 
    Int8       =  256 
    UInt16     =  512 
    Int16      =    4 
    UInt32     =  768 
    Int32      =    8 
    UInt64     = 1280 
    Int64      = 1024 
    Float32    =   16 
    Float64    =   64 
    Float128   = 1536 
    Complex64  =   32 
    Complex128 = 1792 
    Complex256 = 2048 
    RGB24      =  128 

    Map = {
        # Type:  ( NArray , NArrayShort, Pack , Bitpix )
        # Type name = NumPy string = Nifti #define
        Bool:    ( "Bool", "", "c", 1 ),
        Int8:    ( "Int8", "", "b", 8 ),
        UInt8:   ( "UInt8", "", "B", 8 ),
        Int16:   ( "Int16", "i2", "h", 16 ),
        UInt16:  ( "UInt16", "u2", "H", 16 ),
        Int32:   ( "Int32", "i4", "i", 32 ),
        UInt32:  ( "UInt32", "u4", "I", 32 ),
        Float32: ( "Float32", "f4", "f", 32 ),
        Float64: ( "Float64", "f8", "d", 64 ),
        Complex64:  ( "Complex64", "", "ff", 64 ),
        Complex128: ( "Complex128", "", "dd", 128 ),
    }

class RawWriter:

    def __init__(self,filename,gzip=False):
        self.filename = filename
    
        if gzip:
            self.fh = gzip.GzipFile(filename, 'wb')
        else:
            self.fh = file(filename, 'wb')

        # currently expect user to set:
        self.type = NiftiType.Int16
        self.data = None
        self.offset = 0
        self.order = "<"

    def write(self):
        self.fh.close()
        nastr   = self.order + NiftiType.Map[self.type][1]
        dim = self.data.shape
    
        flat = numpy.reshape(self.data, dim[0]*dim[1]*dim[2], order="F")
        mm   = numpy.memmap(self.filename, dtype=nastr, mode="r+",
                offset=self.offset, shape=(dim[0]*dim[1]*dim[2],))
        mm[:] = flat[:]
        del mm

        self.fh = file(self.filename, 'rb+')

        # self.fh.seek(self.offset)
        # packstr = self.order + NiftiType.Map[self.type][2]
        # for e in flat:
        #     bytes = struct.pack(packstr, int(e)) # numpy-related fix
        #     self.fh.write(bytes)
