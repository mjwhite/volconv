#!/usr/bin/env python
#
# Volconv - geometry-aware DICOM-to-NIfTI converter
# DICOM file and series reader
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
import re
import sys
import os
import pydcm
import copy
import orient
import math
import string
import pprint
import fnmatch
from datetime import datetime
import numpy as np
from os.path import basename

ver = map(int, string.split(
    re.sub(r'rc\d+$','',string.split(sys.version)[0]),"."))

if ver < [2,4,0]:
    from sets import Set
    def set(x):
        return Set(x)

def fnmatch_cpt(pathname,pattern):
    pathl = pathname.split(os.path.sep)
    pattl = pattern.split(os.path.sep)

    # cannot possibly match if the pattern is longer than the path
    if len(pattl) > len(pathl):
        return False

    else:
        for i in range(0,len(pattl)):
            if not fnmatch.fnmatch(pathl[i],pattl[i]):
                return False

        return True

def puts(str):
    sys.stderr.write(str)
    sys.stderr.flush()

def roundto(x,sf=3):
    m, e = math.frexp(x)
    m2 = round(m * (10**sf)) / (10**sf)
    return math.ldexp(m2,e)

def lowest(t1,t2):
    """Find consistent lowest of two tuples/lists"""
    compare_len = min(len(t1), len(t2))
    for i in range(0,compare_len):
        if t1[i] < t2[i]:
            return t1
        elif t1[i] > t2[i]:
            return t2

    # if here, identical to compare_len; just pick one
    return t1

def distance(t1,t2):
    """Find vector distance between two tuples"""

    if len(t1) != len(t2):
        raise ValueError("distance: tuples must be the same length")

    r = [0,] * len(t1)
    for i in range(0,len(t1)):
        r[i] = (t1[i] - t2[i]) ** 2
    r = tuple(r)

    return math.sqrt( sum(r) )

def angle(o1,o2):
    """
    Find the angles between two DICOM orientation vectors
    """

    o1 = np.array(o1)
    o2 = np.array(o2)

    o1a = o1[0:3]
    o1b = o1[3:6]
    
    o2a = o2[0:3]
    o2b = o2[3:6]

    norm_a = np.linalg.norm(o1a) * np.linalg.norm(o2a)
    norm_b = np.linalg.norm(o1b) * np.linalg.norm(o2b)

    dot_a = np.dot(o1a,o2a) / norm_a
    dot_b = np.dot(o1b,o2b) / norm_b
    
    if dot_a > 1.0  and dot_a - 1.0 <= np.finfo(dot_a.dtype).eps:
        dot_a = 1.0
    
    if dot_b > 1.0  and dot_b - 1.0 <= np.finfo(dot_b.dtype).eps:
        dot_b = 1.0

    angle_a = np.arccos(dot_a) * (180.0 / np.pi)
    angle_b = np.arccos(dot_b) * (180.0 / np.pi)

    return (angle_a, angle_b)

def plural(n):
    if n == 1:
        return ''
    else:
        return 's'

def unpack1(fmt,str):
    s = struct.calcsize(fmt)
    val = struct.unpack(fmt,str[0:s])
    rest = str[s:]
    return (val, rest)

def fixser(sno):
    snop = re.compile(r'(\d+)(\w*)')
    snom = snop.match(sno)
    snos = "%04d%s" % (int(snom.group(1)), snom.group(2),)
    return snos
    
pat_under = re.compile(r'[\s/^]')
pat_removep = re.compile(r'[^A-Za-z0-9,.;:=%^&()_+-]')

def tidy_protoname(desc):
    """
    Tidy protocol names/descriptions in a uniform way 
    so they can serve as usable filenames
    """

    descc = pat_under.subn('_', desc)[0]
    descc = pat_removep.subn('', descc)[0]
    return descc

def alpha_ser(n):
    letters = [chr(e) for e in range(97,97+25)] # a .. y
    if n <= 24:
        return letters[n]
    elif n >= 25:
        return 'z%04d' % (n,)

def z_ser(n):
    return 'z%04d' % (n,)

class CSA:
    def __init__(self, d):
        self.d = d

class DicomError(ValueError):
    def __init__(self,string,file):
        self.err = string
        self.file = file

class DicomDict:

    __single = {"full":False}

    def __init__(self):
        self.__dict__ = self.__single

        if not self.full:
            self.vr   = {}
            self.mult = {}
            self.name = {}
            self.shortname = {}
            self.byname = {}
            self.readDict()
            self.full = True

    def readDict(self):
        path = os.path.join(pydcm.__path__[0], "dic", "dicomV3.dic")

        dict = file(path, "r")
        sep  = re.compile(r'\s+')
        emp  = re.compile(r'^\s*$')
        for l in dict:
            if emp.match(l): continue
            l    = l.rstrip()
            elts = sep.split(l,4)
            n    = ( int(elts[0],16), int(elts[1],16) )
            self.vr[n]   = elts[2]
            self.mult[n] = elts[3]
            self.name[n] = elts[4]

            shortname = re.sub(r'\(RET\)','',elts[4])
            shortname = shortname.rstrip(' ')
            shortname = shortname.lstrip(' ')
            shortname = re.sub(r'\s+','_',shortname)
            shortname = re.sub(r'[^A-Za-z0-9_]','',shortname)
            shortname = shortname.lower()
            self.shortname[n] = shortname
            self.byname[shortname] = n

        dict.close()
        return self

class DicomDataset:

    def __init__(self):
        self.__dict__['vals'] = {}
        self.__dict__['dict'] = DicomDict()

    def gettag(self,tag):
        if type(tag) == tuple:
            return tag
        elif type(tag) == str:
            try:
                return self.dict.byname[tag]
            except:
                raise AttributeError('Unknown DICOM DE name')

    def __getattr__(self,name):
        try:
            tag = self.gettag(name)
            return self.vals[tag]
        except KeyError:
            raise AttributeError('No value for this DE')
    
    def __getitem__(self,name):
        try:
            tag = self.gettag(name)
            return self.vals[tag]
        except KeyError:
            raise AttributeError('No value for this DE')

    def __setattr__(self,name,value):
        tag = self.gettag(name)
        self.__dict__['vals'][tag] = value
    
    def __setitem__(self,name,value):
        tag = self.gettag(name)
        self.__dict__['vals'][tag] = value

    def __len__(self):
        return len(self.__dict__['vals'])

    def __delattr__(self,name):
        try:
            tag = self.gettag(name)
            vals = self.__dict__['vals']
            del vals[tag]
        except KeyError:
            raise AttributeError('No value for this DE')
    
    def __delitem__(self,name):
        try:
            tag = self.gettag(name)
            vals = self.__dict__['vals']
            del vals[tag]
        except KeyError:
            raise AttributeError('No value for this DE')

class DicomReader:

    def __init__(self, filename, flat=False, sf=5, csa=1, acr=0):
        self.dict = DicomDict()
        self.fn = filename
        self.fh = file(self.fn, "rb")
        self.level = 0
        self.vals = {}
        self.flat = flat
        self.end = "<" # start off little-endian
        self.csa = csa # unused now
        self.acr = acr

        self.csadata = {}

    def checkType(self):
        self.fh.seek(128)
        prefix = self.fh.read(4)
        return (prefix == "DICM")

    def checkACR(self):
        self.fh.seek(0)
        prefix = self.fh.read(2)
        zero = struct.unpack('<H',prefix)[0] # try little-endian

        if zero == 0x0001 or zero == 0x0002 or zero == 0x0003 or zero == 0x0004 or \
           zero == 0x0005 or zero == 0x0006 or zero == 0x0007 or zero == 0x0008:
            return 1

        if zero == 0x0100 or zero == 0x0200 or zero == 0x0300 or zero == 0x0400 or \
           zero == 0x0500 or zero == 0x0600 or zero == 0x0700 or zero == 0x0800:
            self.end = ">" # best guess: switch to big-endian
            return 1

        return False

    def printin(self,str):
        sys.stdout.write (self.level * "    " + str)

    def dump(self, unknown=True, trunc=True):
        self.dumpTree(self.vals, unknown, trunc)

    def getCSA(self, type, key):
        try:
            csa = self.csadata[type]
        except KeyError:

            if type == 'image':
                vf = self.vals[(0x0029,0x1010)]
            elif type == 'series':
                vf = self.vals[(0x0029,0x1020)]

            csa = self.convertCSA2(vf.d)
            self.csadata[type] = csa

        try:
            v = [v["val"] for v in csa[key]["items"] if v["subhdr"][0] > 0]
            return v
        except KeyError:
            return []
        
    def dumpCSA(self,csa,trunc=1):
        for k in csa.keys():
            v = [v["val"] for v in csa[k]["items"] if v["subhdr"][0] > 0]
            value = "[" + (", ".join(v).replace("\n"," ")) + "]"
            if trunc and len(value) > 60 or k == "MrPhoenixProtocol":
                value = value[0:60] + " ..."
            self.printin(": %s = %s\n" % (k, value,))

    def dumpCSAtype(self,type,trunc=1):
        if type == "image":
            csa = self.convertCSA2(self.vals[(0x0029,0x1010)].d)
            self.dumpCSA(csa,trunc)
        elif type == "series":
            csa = self.convertCSA2(self.vals[(0x0029,0x1020)].d)
            self.dumpCSA(csa,trunc)

    def dumpTree(self, tree, unknown=True, trunc=True):
        keys = tree.keys()
        keys.sort()
        for k in keys:
            if k in self.dict.name:
                name = self.dict.name[k]
            else:
                if not unknown: continue
                name = "(Unknown)"

            if isinstance(tree[k],CSA):
                self.printin("%04x|%04x %036s : (CSA)\n" % (k[0], k[1], name))
                self.dumpCSA(self.convertCSA2(tree[k].d))

            elif isinstance(tree[k],dict):
                self.printin("%04x|%04x %036s : (Sequence) =>\n" %
                        (k[0], k[1], name))
                self.level += 1
                self.dumpTree(tree[k], unknown)
                self.level -= 1

            else:
                if trunc:
                    self.printin("%04x|%04x %036s : %-40.40s \n" %
                            (k[0], k[1], name, "%r %s" % (tree[k], type(tree[k]))))
                else:
                    self.printin("%04x|%04x %036s : %s \n" %
                            (k[0], k[1], name, "%r %s" % (tree[k], type(tree[k]))))

    def convertCSA2(self, str):
        csa = {}

        v, str = unpack1("4s", str)
        if v[0] != "SV10":
            return {}
        
        # unused 4 bytes
        v, str = unpack1("4s", str)

        # number of fields
        v, str = unpack1("<I", str)
        n      = v[0]

        # unused 4 bytes (77)
        v, str = unpack1("4s", str)

        for i in range(0,n):

            # field name
            v, str = unpack1("64s", str)
            fnlen = v[0].find("\0")
            name  = v[0][0:fnlen]
            csa[name] = {}

            # vm
            v, str = unpack1("<i", str)
            csa[name]["vm"] = v[0]

            # vr
            v, str = unpack1("4s", str)
            csa[name]["vr"] = v[0].rstrip("\0")

            # syngodt
            v, str = unpack1("<i", str)
            csa[name]["syngodt"] = v[0]
            
            # number of items within field
            v, str = unpack1("<i", str)
            ni = v[0]
            
            # unused 4 bytes (77)
            v, str = unpack1("<i", str)

            csa[name]["nitems"] = ni
            csa[name]["items"] = [{} for i in range(0,ni)]

            for i in range(0,ni):

                # mysterious 16-byte "xx" sub-header
                subhdr, str = unpack1("<4i", str)
                sublen = subhdr[1]
                csa[name]["items"][i]["subhdr"] = subhdr

                # extract value and skip up to next 4-byte multiple
                val   = str[0:sublen]
                vlen  = val.find("\0")
                val   = val[0:vlen].rstrip(" ")
                csa[name]["items"][i]["val"] = val
                nrem  = (4 - (sublen % 4)) % 4
                str   = str[(sublen+nrem):]

        return csa

    def convertVal(self, de, vr, vl, vf):
        if   (vr == "AE"
            or  vr == "AS"
            or  vr == "CS"
            or  vr == "DA"
            or  vr == "DS"
            or  vr == "DT"
            or  vr == "IS"
            or  vr == "LO"
            or  vr == "LT"
            or  vr == "OB"
            or  vr == "OW"
            or  vr == "PN"
            or  vr == "SH"
            or  vr == "ST"
            or  vr == "TM"
            or  vr == "UI"
            or  vr == "UN"
            or  vr == "UT"):

            # vf may be padded with space to even len...
            vf = vf.rstrip()

            # ...or with a NULL
            if len(vf) > 1 and len(vf) % 2 == 0 and vf[-1] == "\x00":
                vf = vf[0:-1]

            # (this de-padding approach is rather simplistic)

            if de in self.dict.mult:
                mult = self.dict.mult[de].split("-")
            else:
                mult = "1"

            def conv(x):
                try:
                    return int(x)
                except:
                    return "unlimited"

            mult = [conv(x) for x in mult]

            if len(mult)==1 and mult[0]==1:
                vf = [vf]
            elif len(mult)==1 and mult[0]>1:
                vf = vf.split("\\", mult[0]-1)
            elif len(mult)>1 and mult[1]!="unlimited":
                vf = vf.split("\\", mult[1]-1)
            else: # unlimited
                vf = vf.split("\\")

            vf = tuple(vf)

        elif  vr == "AT":
            n  = vl / 2
            vf = struct.unpack(self.end+"H"*n, vf)
        elif  vr == "FL":
            n  = vl / 4
            vf = struct.unpack(self.end+"f"*n, vf)
        elif  vr == "FD":
            n  = vl / 8
            vf = struct.unpack(self.end+"d"*n, vf)
        elif  vr == "SL":
            n  = vl / 4
            vf = struct.unpack(self.end+"i"*n, vf)
        elif  vr == "SS":
            n  = vl / 2
            vf = struct.unpack(self.end+"h"*n, vf)
        elif  vr == "UL":
            n  = vl / 4
            vf = struct.unpack(self.end+"I"*n, vf)
        elif  vr == "US":
            n  = vl / 2
            vf = struct.unpack(self.end+"H"*n, vf)
        else:
            fields = "(0x%04x,0x%04x)" % (de[0], de[1])
            raise DicomError("unknown VR %s in %s, giving up on file"%(vr,fields),self.fn)

        if len(vf)==1:
            vf = vf[0]

        return vf

    def readHeader(self):

        if self.checkType():
            self.fh.seek(132)
            implicit=0

        elif self.acr and self.checkACR():
            self.fh.seek(0)
            implicit=1

        else:
            if self.acr: 
                raise DicomError("not a DICOM or (probably) ACR file",self.fn)
            else:
                raise DicomError("not a DICOM file",self.fn)

        try:
            self.vals = self.readFields(implicit=implicit)
        except:
            raise DicomError("failure reading header",self.fn)
        return self

    def readFields(self, maxbytes=0, implicit=0):
        myvals = {}
        startb = self.fh.tell()
        switch_endian = 0
        switch_implicit = 0
        switch_at = 0
        while 1:
            value_start = self.fh.tell()

            if switch_endian and value_start >= switch_at:
                self.end = ">"
                switch_endian = 0
            
            if switch_implicit and value_start >= switch_at:
                implicit = 1
                switch_implicit = 0

            if maxbytes > 0 and value_start >= startb+maxbytes:
                break

            de = self.fh.read(4)   # data element

            # check for EOF
            if de == "":
                break

            # unpack DE and construct hex representation
            de = struct.unpack(self.end+"HH", de)
            dehex = ["%04x" % x for x in de]

            # implicit end of sequence code
            if de == (0xfffe, 0xe0dd):
                self.fh.read(4)    # value length (should be zero)
                break

            # other implicit sequence codes
            if de[0] == 0xfffe:
                self.fh.read(4)    # value length (should be zero)
                continue

            if implicit:
                if de in self.dict.vr:
                    vr = self.dict.vr[de]
                else:
                    vr = 'UN'

                vl = struct.unpack(self.end+"I", self.fh.read(4))[0]
                vltype = "Long"

            else:

                vr = self.fh.read(2)   # value representation
                tmp = struct.unpack(self.end+"H", self.fh.read(2))[0]

                # read length, according to VR
                if vr == "OB" or vr == "OW" or vr == "SQ" or vr == "UN":
                    vl = struct.unpack(self.end+"I", self.fh.read(4))[0]
                    vltype = "Long"
                else:
                    vl = tmp
                    vltype = "Short"

            # simply skip if unknown VR
            if vr == None:
                self.fh.seek(vl, 1)

            else:

                # recurse if sequence VR
                if vr == "SQ" or vl == 0xFFFFFFFF:
                    if vl > 0:
                        self.level += 1
                        vf = self.readFields(maxbytes=vl,implicit=implicit)
                        if self.flat:
                            for de2 in vf.keys():
                                myvals[de2] = vf[de2]
                            vf = "(flattened)"
                        self.level -= 1
                    else:
                        pass

                # otherwise read VF
                else:

                    # for pixel data element, store location and length
                    if de == (0x7fe0, 0x0010):
                        vf = (self.fh.tell(), vl)
                        self.fh.seek(vl, 1)

                    elif de == (0x0029, 0x1010):
                        vf = CSA(self.fh.read(vl))
                    
                    elif de == (0x0029, 0x1020):
                        vf = CSA(self.fh.read(vl))

                    else:
                        vf = self.fh.read(vl)
                        try:
                            vf = self.convertVal(de, vr, vl, vf)
                        except:
                            raise DicomError("VR error, giving up on file %s",self.fn)


                    if de == (0x0002, 0x0000):
                        switch_at = value_start + vf

                    # check transfer syntax
                    if de == (0x0002, 0x0010):
                        if   vf == "1.2.840.10008.1.2": # implicit LE
                            switch_implicit = 1
                        elif   vf == "1.2.840.10008.1.2.1": # explicit LE
                            switch_endian = 0
                        elif vf == "1.2.840.10008.1.2.2": # explicit BE
                            switch_endian = 1
                        else:
                            raise DicomError("unhandled TS %s, giving up on file"%(repr(vf),),self.fn)

                myvals[de] = vf

        # end: while 1
        return myvals

class Entity:
    def __repr__(self):
        return pprint.pformat(self.__dict__)

class DicomSequenceReader:

    def __init__(self, paths, pattern='', flat=False, timehack=False,
            seqinc='', seqexc='', csa=1, acr=0, splitorient=True, single=False, 
            mosaic=None, slice3d=False, sliceinst=False, stackunk=False, sar=False, phase=False,
            fnmatch=None, fnmatch_relative=False, roundorient=True, roundorientthresh=0.2,
            nsubseries=False,
            typeinc='', typeexc=''):
    
        self.files = []
        self.csa = csa
        self.acr = acr
        self.splitorient = splitorient
        self.single = single
        self.mosaic = mosaic
        self.slice3d = slice3d
        self.sliceinst = sliceinst
        self.stackunk = stackunk
        self.sar = sar
        self.phase = phase
        self.roundorient = roundorient
        self.roundorientthresh = roundorientthresh
        self.nsubseries = nsubseries

        self.typeinc = typeinc
        self.typeexc = typeexc

        if pattern == "":
            self.pattern = None
        else:
            self.pattern = re.compile(pattern)
        self.fnmatch = fnmatch
        
        def extendby(basepath,dirname,names):
            if not fnmatch is None:
                if fnmatch_relative:
                    pfnmatch = os.path.join(basepath,self.fnmatch)
                else:
                    pfnmatch = self.fnmatch
            else:
                pfnmatch = None

            for x in names:
                px = os.path.join(dirname,x)
                if not os.path.isdir(px):
                    if (pfnmatch is None) or fnmatch_cpt(px,pfnmatch):
                        if (self.pattern is None) or self.pattern.search(px):
                            self.files.extend([px])

        for path in paths:
            if os.path.isdir(path):
                os.path.walk(path,extendby,path)
            else:
                self.files.append(path)
        
        ## this is neater, but os.walk isn't available in Python 2.2
        #for pn in os.walk(path):
        #    self.files.extend([(pn[0]+os.sep+x) for x in pn[2]])

        self.seqinc = re.compile(seqinc)
        if seqexc != '':
            self.seqexc = re.compile(seqexc)
        else:
            self.seqexc = None

        self.flat = flat
        self.change_name = None
        self.use_exdcm = True
        self.exdcm_path = False
        self.timehack = timehack
        self.show_error_eg = True

    def dumpStudies(self):
        pprint.pprint(self.studies)

    def scanAll(self):
        self.studies  = {}

        n = -1
        repeat = 0
        total = len(self.files)
        self.seriescount = 0
        errors = {}
        visited = {}
        orientations = {}
        warnings = []
        errcount = 0
        ws = re.compile(r' ')
        single_study = None
        single_name = None
        single_ser = None
        single = self.single

        while 1:

            read_header = 0
            
            if not repeat:
                n += 1
                read_header = 1

            if n == len(self.files):
                break

            f = self.files[n]

          # if not self.pattern.search(f): continue

          # if not self.fnmatch is None:
          #     if not fnmatch_cpt(f,self.fnmatch): continue

            repeat = 0
            warnings=[]
            instance_time = 0
            try:
                if read_header:
                    d = DicomReader(f,self.flat,5,self.csa,self.acr).readHeader()

                puts("\rReading: %i/%i (%i warning%s)  "%(n+1,total,errcount,plural(errcount)))

                try:

                    # start with biographical/seq data for filtering
                    try:
                        # study  = d.vals[0x0020,0x0010]
                        study  = d.vals[0x0020,0x000D]
                    except KeyError:
                        study  = "anon"
                    
                    try:
                        name   = d.vals[0x0010,0x0010]
                    except KeyError:
                        name = "anon"

                    if single_study == None:
                        single_study = study
                    
                    if single_name == None:
                        single_name = name

                    if single:
                        name = single_name + "_S"
                        study = single_study + "_S"
                   
                    try:
                        echo  = int(d.vals[0x0018,0x0086])
                    except KeyError:
                        echo  = 1
                    
                    try:
                        te  = float(d.vals[0x0018,0x0081])
                    except (KeyError, ValueError):
                        te  = 0.0

                    try:
                        tr  = float(d.vals[0x0018,0x0080])
                    except (KeyError, ValueError):
                        tr  = 0.0
                    
                    try:
                        flip  = float(d.vals[0x0018,0x1314])
                    except (KeyError, ValueError):
                        flip  = 0.0

                    try:
                        vflip  = d.vals[0x0018,0x1315]
                    except (KeyError, ValueError):
                        vflip  = 'N'

                    try:
                        xdesc   = d.vals[0x0008,0x103e]
                    except KeyError:
                        try:
                            xdesc   = d.vals[0x0018,0x1030]
                        except KeyError:
                            try:
                                xdesc   = d.vals[0x0008,0x1030]
                            except:
                                xdesc   = "unknown"
                
                    # run desc exclusions before we try other parameters (which might fail)
                    if not self.seqinc.search(xdesc):
                        raise DicomError("description didn't match include pattern, skipping file", f)
                    
                    if self.seqexc != None and self.seqexc.search(xdesc):
                        raise DicomError("description matched exclude pattern, skipping file", f)

                    try:
                        image_cmt = d.vals[0x0020,0x4000]
                    except KeyError:
                        image_cmt = None

                    try:
                        patient_cmt = d.vals[0x0010,0x4000]
                    except KeyError:
                        patient_cmt = None

                    # whole image type
                    xtype   = "/".join(d.vals[0x0008,0x0008])

                    # modality-specific type, used to modify the sequence name
                    image_type = d.vals[0x0008,0x0008]
                    if len(image_type) > 2:
                        ximtype = ws.subn('_', image_type[2])[0].lower()
                    else:
                        ximtype = ''

                    # run type exclusions next
                    if self.typeinc != "" and (not self.typeinc.upper() in [e.upper() for e in image_type]):
                        raise DicomError("type didn't match include value, skipping file", f)
                    
                    if self.typeexc != "" and (self.typeexc.upper() in [e.upper() for e in image_type]):
                        raise DicomError("type matched exclude value, skipping file", f)

                    # extract some useful parameters
                    try:
                        ser    = d.vals[0x0020,0x0011]
                    except KeyError:
                        ser    = "0"

                    # one broken GE dataset needed this...
                    ser    = ser.lstrip(' ')

                    if single_ser == None:
                        single_ser = ser

                    if single:
                        ser = single_ser + "S"

                    no_geometry = False
                    if not(d.vals.has_key( (0x0020,0x0032) ) and
                            d.vals.has_key( (0x0020,0x0037) )):
                        if self.stackunk or self.sliceinst:
                            no_geometry = True
                            warnings.append(DicomError("unknown geometry, using naive stacking", f))
                        else:
                            raise DicomError("no geometry, no --stack-unk, skipping file", f)

                    if self.slice3d:
                        pos = [float(x) for x in d.vals[0x0020,0x0032]]
                        orn = [float(x) for x in d.vals[0x0020,0x0037]]
                        k = [
                                orn[1]*orn[5] - orn[2]*orn[4],
                                orn[2]*orn[3] - orn[0]*orn[5],
                                orn[0]*orn[4] - orn[1]*orn[3],
                            ]
                        slice = "%f" % (k[0]*pos[0] + k[1]*pos[1] + k[2]*pos[2],)
                        sliced = d.vals[0x0020,0x0032]
                    elif self.sliceinst or no_geometry:
                        sliced = ["0.0", "0.0", str(float(d.vals[0x0020,0x0013]))]
                        slice = sliced[2]
                    else:
                        try:
                            slice  = d.vals[0x0020,0x1041]  # orthogonal slice location
                            sliced  = d.vals[0x0020,0x0032] # 3D location
                        except KeyError:
                            slice = d.vals[0x0020,0x0032][2] # last elt of 3D slice position (may fail)
                            sliced = d.vals[0x0020,0x0032]   # 3D location

                    try:
                        instance = int(d.vals[0x0020,0x0013])
                    except KeyError:
                        instance = 1
                        warnings.append(DicomError("missing instance number, assuming 1", f))

                    # dtimes = d.vals[0x0020,0x0105]
                    try:
                        time   = d.vals[0x0020,0x0100]
                    except KeyError:
                        try:
                            time = d.vals[0x0020,0x0013] # instance number
                            instance_time = 1
                        except KeyError:
                            time   = "0"

                    rows   = d.vals[0x0028,0x0010]
                    cols   = d.vals[0x0028,0x0011]
                    bytes  = d.vals[0x0028,0x0100] / 8

                    try:
                        res    = \
                            [float(x) for x in (d.vals[0x0028,0x0030] + (d.vals[0x0018,0x0088],))]
                    except KeyError:
                        try:
                            res    = \
                                [float(x) for x in (d.vals[0x0028,0x0030] + (d.vals[0x0018,0x0050],))]
                        except ValueError:
                            res = \
                                [float(x) for x in (d.vals[0x0028,0x0030] + ("1",))]
                            warnings.append(DicomError("unknown slice thickness, assuming 1mm", f))
                        except KeyError:
                            res = [1.0, 1.0, 1.0]
                            warnings.append(DicomError("unknown resolution, assuming 1x1x1mm", f))


                    # find a date - we prefer a study date (fixed for the whole study) over
                    # a series-specific date (unless the study spans midnight...)
                    try:
                        xdate   = d.vals[0x0008,0x0020] # study date
                        if int(xdate) == 0: raise KeyError
                    except KeyError:
                        try:
                            xdate   = d.vals[0x0008,0x0021] # series date
                            if int(xdate) == 0: raise KeyError
                        except KeyError:
                            try:
                                xdate   = d.vals[0x0008,0x0022] # acquisition date
                            except KeyError:
                                xdate = "00000000"

                    # find a series time - we prefer a series time, if possible
                    try:
                        xtime   = d.vals[0x0008,0x0031] # series
                    except KeyError:
                        try:
                            xtime   = d.vals[0x0008,0x0030] # study
                        except KeyError:
                            xtime   = "0000"

                    try:
                        study_time   = d.vals[0x0008,0x0030] # study
                    except:
                        study_time = "0000"

                    try:
                        study_date   = d.vals[0x0008,0x0020] # study date
                        if int(study_date) == 0: raise KeyError
                    except:
                        study_date = "00000000"

                    # store the acquisition time (seems to vary hugely: on some scanners it's
                    # the series time, some it's the volume time, and some it's the actual
                    # slice time within a multi-slice sequence)
                    try:
                        dtime   = d.vals[0x0008,0x0032] # acquisition
                    except KeyError:
                        dtime   = xtime # copy series/study time

                    if self.sliceinst or no_geometry:
                        orientt  = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0)
                    else:
                        orientt  = tuple([float(x) for x in d.vals[0x0020,0x0037]])

                    try:
                        intercept = float(d.vals[0x0028,0x1052])
                        slope     = float(d.vals[0x0028,0x1053])
                    except KeyError:
                        intercept = 0.0
                        slope     = 1.0

                    # table position
                    try:
                        table = [int(e) for e in d.vals[0x0019,0x1014].split('\\')]
                    except KeyError:
                        table = None

                    # actual pixel data
                    pixels = d.vals[0x7fe0,0x0010]

                    mosaicid = None

                    # forced mosaic size
                    if self.mosaic:
                        mosaic = [self.mosaic]

                    # definative check (forces CSA "image" parse, which is slow)
                    elif self.csa:
                        mosaic = d.getCSA("image", "NumberOfImagesInMosaic")
                    
                    # attempt fast check
                    else:
                        try:
                            s = d.vals[0x0008,0x0008]
                            if ' '.join(s).upper().find("MOSAIC") != -1:
                                mosaic = d.getCSA("image", "NumberOfImagesInMosaic")
                            else:
                                mosaic = []
                        except KeyError:
                            mosaic = []

                    # StartFMRI "images" appear to have all the CSA headers of the real thing, but with
                    # a small dummy image; so you don't _really_ want to try and unpack them as mosaics
                    # Check for this by looking for "DUMMY IMAGE" in the type string; this doesn't
                    if mosaic != []:
                        s = d.vals[0x0008,0x0008]
                        if ' '.join(s).upper().find("DUMMY IMAGE") != -1:
                            mosaic = []
                            warnings.append(
                                    DicomError("not unpacking mosaic for dummy image", f))
                    
                    # Siemens diffusion direction; this will be a function of time
                    # for Siemens at least, this is in the DICOM patient co-ordinate system
                    try:
                        if self.csa:
                            diff = [float(x) for x in d.getCSA("image","DiffusionGradientDirection")]
                        else:
                            diff = d.vals[0x0019,0x100e]

                            # this should be a float vector, but when LEI images are sent PACS->3T->physics
                            # group 19 becomes entirely "UN" VRs, so it's garbage.  Typically we'll have
                            # the CSA header already by this point for mosaics, so use that instead.
                            if type(diff) == str:
                                try:
                                    diff = [float(x) for x in d.getCSA("image","DiffusionGradientDirection")]
                                    warnings.append(
                                            DicomError("diffusion vector has bad type, using CSA instead", f))
                                except:
                                    diff = []
                                    warnings.append(
                                            DicomError("diffusion vector has bad type, no CSA available", f))
                    except KeyError:
                        diff = []

                    # SAR
                    if self.sar:
                        sar_values        = d.getCSA("series","SARMostCriticalAspect")
                        sar_body_pred     = d.vals[0x0018,0x1316]
                        sar_most_crit     = d.getCSA("series","RFSWDMostCriticalAspect")
                        sar_mode          = d.getCSA("series","RFSWDOperationMode")
                        sar = {
                                "values": [float(x) for x in sar_values],
                                "body":   float(sar_body_pred),
                                "most_crit": sar_most_crit[0],
                                "mode": int(sar_mode[0]),
                                }
                    else:
                        sar = None

                    # Phase encode direction
                    if self.phase:
                        phase_direction     = d.vals[0x0018,0x1312]
                        phase_positive      = int(d.getCSA("image","PhaseEncodingDirectionPositive")[0])

                        if phase_direction == "ROW":
                            phase_axis = "i"
                        elif phase_direction == "COL":
                            phase_axis = "j"

                        if phase_positive == 0:
                            phase_axis = "-" + phase_axis

                        phase = {
                                "direction":    phase_direction,
                                "positive":     phase_positive,
                                "axis":         phase_axis,
                                }
                    else:
                        phase = None

                    # B value
                    try:
                        if self.csa:
                            bval = float(d.getCSA("image","B_value")[0])
                        else:
                            bval = float(d.vals[0x0019,0x100c])
                    except KeyError:
                        bval = None

                    if mosaic != []:
                        nmos = int(mosaic[0])

                        # get visit count
                        try:
                            seen = visited[f]
                        except KeyError:

                            # warn the user, this is guessed-from-data code
                            warnings.append(
                                DicomError("mosiac is not standards-based, beware geometry", f))

                            visited[f] = 0
                            seen = 0

                        # do we need to visit again?
                        if (seen+1) < nmos:
                            repeat = 1

                        # update visit count
                        visited[f] = visited[f] + 1
                        
                        # actual image matrix
                        fac = math.ceil(math.sqrt(nmos))
                        mrows = rows
                        mcols = cols
                        rows = mrows / fac
                        cols = mcols / fac

                        # row and column co-ordinates of this slice
                        rpos = int(seen) / int(fac)
                        cpos = int(seen) % int(fac)

                        # find normal vector
                        k = [0.0, 0.0, 0.0]
                        i = orientt[0:3]
                        j = orientt[3:]
                        k[0] = i[1]*j[2] - i[2]*j[1]
                        k[1] = i[2]*j[0] - i[0]*j[2]
                        k[2] = i[0]*j[1] - i[1]*j[0]

                        # calculate position of actual first slice (not top of mosaic)
                        truepos = [0.0, 0.0, 0.0]
                        colcor = (float(mcols) - float(cols))/2.0
                        rowcor = (float(mrows) - float(rows))/2.0
                        truepos[0] = float(sliced[0]) + i[0]*res[0]*colcor + j[0]*res[1]*rowcor
                        truepos[1] = float(sliced[1]) + i[1]*res[0]*colcor + j[1]*res[1]*rowcor
                        truepos[2] = float(sliced[2]) + i[2]*res[0]*colcor + j[2]*res[1]*rowcor

                        # calculate position of *this* slice
                        fseen = float(seen)
                        spacing = float(d.vals[0x0018,0x0088])
                        slice = str(float(slice) + spacing * fseen)
                        sliced = (
                            truepos[0] + k[0] * spacing * fseen,
                            truepos[1] + k[1] * spacing * fseen,
                            truepos[2] + k[2] * spacing * fseen,
                        )

                        # store information for unpacking mosaic element
                        mosaicid = Entity()
                        mosaicid.mrows = mrows
                        mosaicid.mcols = mcols
                        mosaicid.n = seen
                        mosaicid.rpos = rpos
                        mosaicid.cpos = cpos

                    # Mosaic images:
                    # - calculate and store image region
                    # - update rows, columns, and image position (how?)
                    # - put file back on the list if needed
                    # - possibly do some orientation-related magic (what?)

                except KeyError, e:
                    fields = "(0x%04x,0x%04x)" % (e[0][0], e[0][1])
                    raise DicomError("missing element %s, skipping file"%(fields,), f)

                if not (study, name) in self.studies:
                    self.studies[study, name] = {}

                v = self.studies[study, name]
               
                if not (study, name) in orientations.keys():
                    orientations[study,name] = {}
                if not ser in orientations[study,name].keys():
                    orientations[study,name][ser] = {}

                if self.roundorient:
                    close_enough = None
                    n_close_enough = 0
                    n_nearly_close_enough = 0

                    # Find any known orientation that's close enough to the current one
                    #
                    # -- note that the behaviour of this strategy when several
                    # slices have nearly-close-enough orientations becomes
                    # potentially dependent on the order the files are read, in
                    # two ways.  When a second just-about-distinct known
                    # orientation is added, other slices may fall within the
                    # rounding error of both, so their assignment can depend on
                    # whether they are read before or after it.  And as further
                    # slices are assigned, the lower_exact value may change,
                    # which might affect which other slices fall within range.
                    #
                    # -- in this edge case, we can't know a priori whether the
                    # user really wanted the volumes split or not, so if this
                    # happens, the user will just need to adjust
                    # --round-threshold up or down to get what they wanted.
                    #
                    # -- but to help the user notice that it's happened, we
                    # check for two warning flags: orientations marked as
                    # different but within double the threshold, and slices
                    # which could be assigned to more than one volume.
                    #
                    # -- we never want to merge orientations for which the image
                    # types are different (messes up GE vol+projection series)

                    for known, known_imtype in orientations[study,name][ser].keys():
                        a1, a2 = angle(orientt,known)

                        if known_imtype != ximtype:
                            continue

                        if a1 < self.roundorientthresh and a2 < self.roundorientthresh:
                            close_enough = known
                            n_close_enough += 1

                        if a1 < self.roundorientthresh*2.0 and a2 < self.roundorientthresh*2.0:
                            n_nearly_close_enough += 1

                        if n_nearly_close_enough > n_close_enough:
                            warnings.append(
                                    DicomError("orientation merge had near miss (< 2*threshold)", f))

                    if close_enough is not None:
                        old_exact = close_enough
                        new_exact = orientt
                        
                        if n_close_enough > 1:
                            warnings.append(
                                    DicomError("orientation merge slice assignment is ambiguous", f))

                        if old_exact != new_exact:
                            lower_exact = lowest(old_exact, new_exact)

                            orientt = lower_exact

                            suff = orientations[study,name][ser][old_exact,ximtype]
                            del orientations[study,name][ser][old_exact,ximtype]
                            orientations[study,name][ser][lower_exact,ximtype] = suff

                def suffix(n):
                    if n == 0:
                        return ''
                    else:
                        return 'o%d' % (n,)

                if self.splitorient:
                    if no_geometry:
                        sersuff = 'unk'
                    elif orientations[study,name][ser].has_key((orientt,ximtype)):
                        sersuff = orientations[study,name][ser][orientt,ximtype]
                    else:
                        sersuff = suffix(len(orientations[study,name][ser]))

                    orientations[study,name][ser][orientt,ximtype] = sersuff
                    ser = ser + sersuff

                # updates to v[ser] must happen *after* sersuff has been added
                if self.roundorient and \
                   close_enough is not None and \
                   old_exact != new_exact:

                        del v[ser].orient[old_exact]
                        v[ser].orient[lower_exact] = True

                # SPM8 and SPM12 write special desc fields; 
                # look for "descrip = " in the SPM source. For example:
                # "3T 3D RM TR=22.5ms/TE=11.2ms/FA=20deg/SO=no 01-Dec-2012 12:01:01.123"

                try:
                    scanoptions = str(d.vals[0x0018,0x0022])
                    if scanoptions == "":
                        scanoptions = "no"
                except KeyError:
                    scanoptions = "no"

                if "MOSAIC" in d.vals[0x0008,0x0008]:
                    mosaic = " Mosaic"
                else:
                    mosaic = ""

                adate = d.vals[0x0008,0x0022] # AcquisitionTime
                atime = d.vals[0x0008,0x0032] # AcquisitionDate

                # AcquisitionTime can be different for different slices
                # (eg with a 2D acquisition) so must be store per s/t/e

                dt = datetime.strptime(adate,"%Y%m%d")

                m = re.match(r'(\d{2})(\d{2})([0-9\.]+)$',atime)
                time_hr = m.group(1)
                time_mn = m.group(2)
                time_sc = m.group(3)

                try:
                    descrip = '%gT %s %s TR=%gms/TE=%gms/FA=%gdeg/SO=%s %s %s:%s:%.5g%s' % \
                    (
                        float(d.vals[0x0018,0x0087]), # MagneticFieldStrength
                        str(d.vals[0x0018,0x0023]), # MRAcquisitionType
                        re.sub(r'\s+','', str(d.vals[0x0018,0x0020])), # ScanningSequence
                        float(tr),
                        float(te),
                        float(flip),
                        scanoptions,
                        dt.strftime("%m-%b-%Y"),
                        time_hr,
                        time_mn,
                        float(time_sc),
                        mosaic,
                    )
                except KeyError:
                    descrip = "missing"
                
                if not ser in v:
                    v[ser] = Entity()
                    v[ser].echoes = {}
                    v[ser].te     = {}
                    v[ser].slices = {}
                    v[ser].slicesd = {}
                    v[ser].times  = {} # normally from instance numbers
                    v[ser].dtimes = {} # dynamic time (per-ser, per-vol, or per-slice)
                    v[ser].file   = {}
                    v[ser].end    = {}
                    v[ser].pixels = {}
                    v[ser].rescale = {}
                    v[ser].mosaic = {}
                    v[ser].diff   = {}
                    v[ser].bval   = {}
                    v[ser].descrip = {}
                    v[ser].shape  = (cols, rows)
                    v[ser].res    = res
                    v[ser].desc   = xdesc
                    v[ser].type   = xtype
                    v[ser].date   = xdate # study/series date
                    v[ser].time   = xtime # series/study time
                    v[ser].stdate = study_date
                    v[ser].sttime = study_time
                    v[ser].sar    = sar
                    v[ser].phase  = phase
                    v[ser].imtype = ximtype
                    v[ser].tr     = tr
                    v[ser].flip   = flip
                    v[ser].vflip  = vflip
                    v[ser].table  = table
                    v[ser].instance_time = instance_time
                    v[ser].patient_cmt = patient_cmt
                    v[ser].image_cmt = image_cmt
                    v[ser].instance = instance
                    self.seriescount += 1
                    v[ser].orient = {}

                # this allows us to catenate multiple orientations
                # (in which case we will ignore orientation data!)
                # (these will, of course, be meaningless if resliced...)
                v[ser].orient[orientt] = True

                # record the smallest instance number for each (sub-)series
                # (this is just used as a sort key)
                if v[ser].instance > instance:
                    v[ser].instance = instance

                sliceoff = 10000.0 * (len(v[ser].orient)-1)
                sliceind = sliceoff + float(slice)

                v[ser].slices[sliceind] = True
                v[ser].slicesd[sliceind] = sliced
                v[ser].echoes[echo] = True
                v[ser].te[echo] = te
                v[ser].times[time]   = True
                v[ser].file[sliceind,time,echo] = f
                v[ser].end[sliceind,time,echo] = d.end
                v[ser].pixels[sliceind,time,echo] = pixels
                v[ser].rescale[sliceind,time,echo] = (intercept, slope)
                v[ser].mosaic[sliceind,time,echo] = mosaicid
                v[ser].dtimes[sliceind,time,echo] = dtime
                v[ser].descrip[sliceind,time,echo] = descrip
                v[ser].diff[time] = diff
                v[ser].bval[time] = bval
                
                for w in warnings:
                    errcount += 1
                    if not errors.has_key(w.err):
                        errors[w.err] = [w.file]
                    else:
                        errors[w.err].append(w.file)

            except DicomError, d:
                errcount += 1
                if not errors.has_key(d.err):
                    errors[d.err] = [d.file]
                else:
                    errors[d.err].append(d.file)
                continue


        # Rename orientation sub-series where possible
        #
        # -- if the series has only one orientation block, no suffix
        # will be added (most common case)
        #
        # -- if all orientations can be given meaningful names (axi, sag,
        # cor, obl), do this (handles localizers and most volumes with a
        # single transverse reconstruction sensibly, which is the next
        # most common case)
        # 
        # -- if not, and there are <= 25 orientations in total, they
        # will be suffixed [a-y] in instance number order, including the
        # baseline series (this is a good choice to handle multi-slab
        # sequences in the spine, for example)
        #
        # -- if there are more than 26 orientations in total, or the
        # user specified --n-subseries, all will be labelled zNNNN (this is
        # the right general solution for tumbling-CoW, MRV, etc)
        #
        # FUTURE: the orientation splitting mechanism should be generalized to
        # handle multiple slabs which *are* parallel
        #
        if self.splitorient:
            for studyk in orientations.keys():
                study = orientations[studyk]
                v     = self.studies[studyk]
                for serk in study.keys():
                    ser = study[serk]
                    if len(ser) > 1:

                        subseries = []
                        for orientt, ximtype in ser:
                            suff = ser[orientt, ximtype]
                            subser = serk + suff
                            serdata = v[subser]

                            subseries.append(
                                    [serdata.instance, orientt, ximtype, serk, suff, subser])

                        def loc_cmp(a,b):
                            return cmp(a[0], b[0])
                        subseries.sort(loc_cmp)

                        mapping = {}

                        if self.nsubseries:
                            short_names = False
                        else:
                            used = {}
                            short_names = True
                            for inst, orientt, ximtype, serk, suff, subser in subseries:
                                orient_str = \
                                    orient.OrientedImage(None,serdata.res,[orientt]).findOrient('short')
                                if not used.has_key(orient_str):
                                    mapping[subser] = serk + orient_str
                                    used[orient_str] = True
                                else:
                                    short_names = False

                        if not short_names:
                            if self.nsubseries:
                                make_ser = z_ser
                            else:
                                make_ser = alpha_ser
                        
                            mapping = {}

                            for i in range(0,len(subseries)):
                                inst, orientt, ximtype, serk, suff, subser = subseries[i]
                                mapping[subser] = serk + make_ser(i)

                        for k in mapping.keys():
                            v[mapping[k]] = v[k]
                            del v[k]

        # If no temporal position identifier was available, instance
        # numbers are recorded instead: so it looks like each slice is
        # at a separate time.  Here will collapse these back into
        # volumes, by setting the time on each block of pixel data to
        # the instance number normalized to the number of time-points
        # it looks like we *should* have.
        # 
        # This will happen routinely for non-dynamic series, so we
        # only want to raise a warning where actual varying dynamic
        # times will be generated.
        #
        for s in self.studies.values():
            for e in s.values():

                if e.instance_time:
                    new_file = {}
                    new_end = {}
                    new_pixels = {}
                    new_times = {}
                    new_dtimes = {}
                    new_rescale = {}
                    new_mosaic = {}
                    new_diff = {}
                    new_bval = {}
                    new_descrip = {}

                    stes = e.file.keys()
                    slices = len(e.slices)
                    echoes = len(e.echoes)
                    numerictimes = [int(x) for x in e.times.keys()]
                    numerictimes.sort()
                    times = len(numerictimes)

                    # accumulate a set of distinct (s,e) groups
                    times_se = {}
                    for ste in stes:
                        slice, time, echo = ste
                        if not times_se.has_key((slice,echo)):
                            times_se[slice,echo] = [int(time)]
                        else:
                            times_se[slice,echo].append(int(time))

                    # how many time-points?
                    nt_set = set([len(times_se[se]) for se in times_se.keys()])
                    nt     = max(nt_set)

                    # this means missing planes - not all se values have same number of t values;
                    # with instance-order stacking, the missing planes will all migrate to the
                    # later t values, which may mess up ordering, too
                    if len(nt_set) > 1:
                        warnings.append(
                            DicomError("missing planes in instance order, gaps may be assigned to wrong volume",
                            e.file[e.file.keys()[0]]))

                    # sort the lists of time-points, build a mapping table, find their intervals
                    # (mapping assumes that, if an se combination appears multiple times, their
                    # volume membership is just the order of their instance numbers; this will work
                    # whether adjacent images are the same slice or the same time-point)
                    times_map = {}
                    all_deltas = []
                    deltas_se = {}
                    for se in times_se.keys():
                        times_se[se].sort()
                        deltas = []
                        t = times_se[se]

                        if len(times_se[se]) > 1:
                            for i in range(0,len(times_se[se])-1):
                                deltas.append( t[i+1] - t[i] )
                                deltas_se[se] = deltas
                        all_deltas.extend(deltas)

                        for ni in range(0,len(times_se[se])):
                            times_map[times_se[se][ni]] = ni

                    # find the modal delta = length of a (s,e) group in instance numbers,
                    # to check for consistency of ordering and/or missing planes
                    all_deltas.sort()
                    if len(all_deltas) == 0:
                        assert(nt == 1)
                        dt = None
                    else:
                        hist = [(i, all_deltas.count(i)) for i in set(all_deltas)]
                        hist.sort(lambda a,b: cmp(b[1],a[1]))
                        dt = hist[0][0]

                    if len(set([tuple(elt) for elt in deltas_se.values()])) > 1:
                        warnings.append(
                            DicomError("instance spacing inconsistent, multi-volume slice assignment may be wrong",
                                e.file[sorted(e.file.keys())[0]]))
                    elif len(set(all_deltas)) > 1:
                        warnings.append(
                            DicomError("instance spacing not constant, series probably has multiple volume axes",
                                e.file[sorted(e.file.keys())[0]]))

                    if nt != times:
                    
                        for ste in stes:
                            slice, time, echo = ste

                            # map instance numbers to actual volume times
                            vol_time = str( times_map[int(time)] )

                            new_file[slice, vol_time, echo] = e.file[ste]
                            new_end[slice, vol_time, echo] = e.end[ste]
                            new_pixels[slice, vol_time, echo] = e.pixels[ste]
                            new_rescale[slice, vol_time, echo] = e.rescale[ste]
                            new_mosaic[slice, vol_time, echo] = e.mosaic[ste]
                            new_dtimes[slice, vol_time, echo] = e.dtimes[ste]
                            new_descrip[slice, vol_time, echo] = e.descrip[ste]
                            new_diff[vol_time] = e.diff[time]
                            new_bval[vol_time] = e.bval[time]

                            new_times[vol_time] = True
                    
                            # if we're generating multiple time-points, warn the user on this file
                            if nt > 1:
                                warnings.append(DicomError("guessing times from instance numbers", 
                                    e.file[ste]))
                        
                        e.times = new_times
                        e.file = new_file
                        e.end = new_end
                        e.pixels = new_pixels
                        e.rescale = new_rescale
                        e.mosaic = new_mosaic
                        e.dtimes = new_dtimes
                        e.descrip = new_descrip
                        e.diff = new_diff
                        e.bval = new_bval


        # work out which volumes are missing slices
        for s in self.studies.values():
            for e in s.values():
                slices = len(e.slices)
                echoes = len(e.echoes)
                times  = len(e.times)
                stes = e.file.keys()
                sl_count = {}

                for ste in stes:
                    slice, time, echo = ste
                    if sl_count.has_key((time,echo)):
                        sl_count[time,echo] += 1
                    else:
                        sl_count[time,echo] = 1
                
                e.missing = {}
                for te in sl_count.keys():
                    missed = slices - sl_count[te]
                    e.missing[te] = missed
                
                if len([i for i in e.missing.values() if i>0]) > 0:
                    warnings.append(DicomError("missing slices in volumes generated from series", 
                        e.file[ste]))

        # print actual warnings
        for w in warnings:
            errcount += 1
            if not errors.has_key(w.err):
                errors[w.err] = [w.file]
            else:
                errors[w.err].append(w.file)

        puts("\rRead: %i/%i (%i warning%s)     \n"%(n,total,errcount,plural(errcount)))
       
        for k in errors.keys():
            puts("Warning: %s (repeated %d time%s)\n"%(k,len(errors[k]),plural(len(errors[k]))))
            if self.show_error_eg: puts("     eg: %s\n"%(errors[k][0],))

    def volumecount(self):
        volumes = 0
        for study in self.studies:
            for series in self.studies[study]:
                volumes += len(self.studies[study][series].times) * \
                    len(self.studies[study][series].echoes)
        return volumes

    def dump(self,alias=None):
        for study in self.studies:
            print "Study:", study
            ksort = self.studies[study].keys()[:]
            ksort.sort(lambda a,b: cmp(fixser(a),fixser(b)))
            for k in ksort:
                e = self.studies[study][k]

                missing_list = e.missing.values()
                if len([i for i in missing_list if i>0]) > 0:
                    missing = " (%d/%d gaps)" %  \
                        (sum(missing_list), len([i for i in missing_list if i > 0]))
                else:
                    missing = ""

                aliasname = ""
                if alias:
                    aliasname = alias.match(study[0],study[1],k)
                    if aliasname is None:
                        aliasname = ""
                    else:
                        if aliasname[1] is None:
                            aliasname = aliasname[0]
                        else:
                            aliasname = aliasname[0] + "-" + str(aliasname[1])
                        aliasname = " {" + aliasname + "}"
                
                description = re.compile(r' ').subn('_', e.desc)[0]
                orient_str = orient.OrientedImage(None,e.res,e.orient.keys()).findOrient()
                print "  Series:", k, "->", \
                        ("(%d x %d x %d)" % (e.shape[0], e.shape[1], len(e.slices),)), \
                        len(e.times), len(e.echoes), orient_str, "[" + description + "]", \
                        e.imtype + missing + aliasname
                
            # print "  Series:", k, "->", e.shape, len(e.slices), "[", \
            #     min(e.slices), max(e.slices), "]", len(e.times), "[", \
            #     min(e.times), max(e.times), "]"

    def anonymize(self, newname="Anonymous", use_exdcm=False):
        self.change_name = newname
        self.use_exdcm = use_exdcm
            
    def exdcmpath(self, exdcmpath):
        self.exdcm_path = exdcmpath

    def interval(self, e):
        if len(e.times) == 1:
            return 0.0
        
        tmptimes = [(float(x),x) for x in e.times.keys()]
        tmptimes.sort(lambda a,b: cmp(a[0],b[0]))

        # indices of first two time-points
        n0 = 0
        n1 = 1
        time0 = tmptimes[n0][1]
        time1 = tmptimes[n1][1]
        nt    = len(tmptimes)

        # chose a random slice and echo (doesn't matter which)
        aslice = e.slices.keys()[0]
        aecho  = e.echoes.keys()[0]

        while 1:
            try:
                interval = float(e.dtimes[aslice, time1, aecho]) - \
                        float(e.dtimes[aslice,time0,aecho])
                break
            except:
                n0 += 1
                n1 += 1

                if n1 == nt:
                    interval = 0.0
                    break

                time0 = tmptimes[n0][1]
                time1 = tmptimes[n1][1]

        return interval

    def toxml(self,filenames={}):
        """
        XML output has been deprecated for a while, and no longer contains as much
        information as the JSON header below.
        """
        out = ""
        for study in self.studies:
            out += "<study id=\"%s\">\n" % study[0]
            out += "<id>%s</id>\n" % study[0]
            if self.change_name:
                out += "<name>%s</name>\n" % self.change_name
            else:
                out += "<name>%s</name>\n" % study[1]
            ksort = self.studies[study].keys()[:]
            ksort.sort(lambda a,b: cmp(fixser(a),fixser(b)))
            for k in ksort:
                e = self.studies[study][k]
                out += "<series id=\"%s\">\n" % k
                out += "  <id>%s</id>\n" % k
                out += "  <rows>%d</rows>\n" % e.shape[0]
                out += "  <cols>%d</cols>\n" % e.shape[1]
                out += "  <slices>%d</slices>\n" % len(e.slices)
                out += "  <times>%d</times>\n" % len(e.times)
                out += "  <echoes>%d</echoes>\n" % len(e.echoes)
                out += "  <flip var=\"%s\">%g</flip>\n" % (e.vflip, e.flip)
                out += "  <reptimes>%g</reptimes>\n" % e.tr

                # if len(e.echoes) > 1:
                tmpechoes = [("%g"% (float(e.te[x]),)) for x in e.echoes.keys()]
                out += "  <echotimes>%s</echotimes>\n" % " ".join(tmpechoes)
                
                if len(e.times) > 1:
                    out += "  <interval>%.4g</interval>\n" % (self.interval(e),)

                out += "  <desc>%s</desc>\n" % e.desc
                out += "  <type>%s</type>\n" % e.type
                if self.use_exdcm:
                    if self.exdcm_path:
                        fn = e.file[sorted(e.file.keys())[0]]
                    else:
                        fn = basename(e.file[sorted(e.file.keys())[0]])
                    out += "  <exdcm>%s</exdcm>\n" % fn
                if filenames.has_key((study[0],study[1],k)):
                    fn = basename(filenames[study[0],study[1],k])
                    if fn.endswith('.nii'):
                        out += "  <nii>%s</nii>\n" % fn
                    elif fn.endswith('.gipl'):
                        out += "  <gipl>%s</gipl>\n" % fn
                out += "  <date>%s</date>\n" % e.date
                out += "  <time>%s</time>\n" % e.time
                out += "</series>\n"
            out += "</study>\n"
        return(out)
    
    def tojson(self,filenames={},alias=None,axes={}):
        """
        The convention for JSON output from volconv is that all geometric information
        is in the DICOM space: [x y z] in the DICOM LPS co-ordinate system, and any
        grid-relative directions like diffusiongrid are in the [i j normk] system,
        regardless of actual slice order, rotation etc.

        The "axes" field has the mapping between DICOM [i j normk] and Nifti [I J K]
        voxel indices.  This depends on the conversion (stack, flip, re-orient) parameters.

        The world mapping is always the same anyway: [X Y Z] = [-x -y z]

        Fields with the suffix _out are mapped into the output spaces [I J K] or [X Y Z].
        It might be even nicer to suffix _xyz, _ijk, _XYZ or _IJK to all geometric
        data in a future release, but at the moment that would break compatibility.
        """

        out = "["
        firstser = True
        for study in self.studies:
            if not firstser:
                out += ',\n'
            else:
                firstser = False
            out += '{\n'
            out += '  "type": "study",\n'
            out += '  "id": "%s",\n' % study[0]
            if self.change_name:
                out += '  "name": "%s",\n' % self.change_name
            else:
                out += '  "name": "%s",\n' % study[1]

            out += '  "series": [\n'
            ksort = self.studies[study].keys()[:]
            ksort.sort(lambda a,b: cmp(fixser(a),fixser(b)))
            first = True
            for k in ksort:

                if alias:
                    if not alias.match(study[0], study[1], k):
                        continue

                e = self.studies[study][k]
                if not first:
                    out += ',\n'
                else:
                    first = False

                out += '  {\n'
                out += '    "type": "series",\n'
                out += '    "id": "%s",\n' % k
                out += '    "rows": %d,\n' % e.shape[0]
                out += '    "cols": %d,\n' % e.shape[1]
                out += '    "slices": %d,\n' % len(e.slices)
                out += '    "times": %d,\n' % len(e.times)
                out += '    "echoes": %d,\n' % len(e.echoes)
                out += '    "flip_var": "%s",\n' % e.vflip
                out += '    "flip": %g,\n' % e.flip
                out += '    "reptimes": [%g],\n' % e.tr
                
                tmpechoes = [("%g"% (float(e.te[x]),)) for x in e.echoes.keys()]
                out += '    "echotimes": [%s],\n' % (', '.join(tmpechoes))

                if not e.table is None:
                    out += '    "table": [%s],\n' % (', '.join([str(x) for x in e.table]),)

                if not e.patient_cmt is None:
                    out += '    "patient_cmt": "%s",\n' % e.patient_cmt
                
                if not e.image_cmt is None:
                    out += '    "image_cmt": "%s",\n' % e.image_cmt

                if not e.sar is None:
                    out += '    "sar": {\n'
                    out += '        "mode": %d,\n' % (e.sar['mode'],)
                    out += '        "most_crit": "%s",\n' % (e.sar['most_crit'],)
                    out += '        "value_lim":  %g,\n' % (e.sar['values'][0],)
                    out += '        "value_1":  %g,\n' % (e.sar['values'][1],)
                    out += '        "value_2":  %g,\n' % (e.sar['values'][2],)
                    out += '        "value_body": %g\n' % (e.sar['body'],)
                    out += '    },\n'

                if not e.phase is None:
                    out += '    "phase": {\n'
                    out += '        "axis": "%s",\n' % (e.phase['axis'],)

                    if axes.has_key((study[0],study[1],k)):
                        a = axes[study[0],study[1],k]
                        out_axis = orient.map_axis(e.phase['axis'],a)
                        out += '        "axis_out": "%s",\n' % (out_axis,)
                    
                    out += '        "direction": "%s",\n' % (e.phase['direction'],)
                    out += '        "positive": %d\n' % (e.phase['positive'],)

                    out += '    },\n'
                
                if len(e.times) > 1:
                    out += '    "interval": %.4g,\n' % (self.interval(e),)
                
                if e.bval.values()[0] != None:
                    t = [int(x) for x in e.times]; t.sort(); t = [str(x) for x in t]

                    # diffusion in the DICOM [x,y,z] co-ordinates
                    out += '    "diffusion": [\n'
                    first=1
                    for tn in t:
                        if first==0:
                            out += ',\n'
                        first=0
                        if len(e.diff[tn]) < 3:
                            out += '        [%g, null]' % (e.bval[tn])
                        else:
                            out += '        [%g, [%g, %g, %g]]' % (e.bval[tn], e.diff[tn][0], e.diff[tn][1], e.diff[tn][2])
                    out += '\n    ],\n'
                    
                    # diffusion aligned to the DICOM [i,j,k] image grid
                    #
                    # (NB we use OrientedImage here re-initialized from the original DICOM
                    # orientation field for this series)
                    out += '    "diffusiongrid": [\n'
                    o = orient.OrientedImage(None,e.res,e.orient.keys())
                    first=1
                    for tn in t:
                        if first==0:
                            out += ',\n'
                        first=0
                        if len(e.diff[tn]) < 3:
                            out += '        [%g, null]' % (e.bval[tn])
                        else:
                            egrid = o.dcm_to_grid(e.diff[tn])
                            out += '        [%g, [%g, %g, %g]]' % (e.bval[tn], egrid[0], egrid[1], egrid[2])
                    out += '\n    ],\n'

                out += '    "desc": "%s",\n' % e.desc
                out += '    "type": "%s",\n' % e.type

                # Axes mappings:
                #
                # This records the DICOM-to-output axis mappings for the conversion.  We only allow
                # for very simple mappings, flips and 90-degree rotations: a major point of volconv
                # is to preserve the image and DICOM co-ordinate systems.
                #
                # The [i j k] system is the image array: [0 0 0] is the corner voxel of the lowest
                # slice.  [I J K] is the output image array with [0 0 0] the first voxel stored.
                # [i j k] and [I J K] indices are always positive; if an axis is flipped, the corner
                # voxel moves to the opposite end of that axis.
                #
                # The [x y z] system is the patient co-ordinate system: [0 0 0] is the DICOM origin.
                # [X Y Z] of the output format will be an appropriate patient/world co-ordinate
                # system with the same origin and orthogonal axes.  The anatomical meaning of DICOM
                # axes is preserved if the output format defines them (eg DICOM LPS -> Nifti RAS).

                if axes.has_key((study[0],study[1],k)):
                    a = axes[study[0],study[1],k]

                    # image axis mappings: depends on reorient/flip parameters
                    out += '    "grid_axes_map": ["%s", "%s", "%s"],\n' % \
                        (a[0], a[1], a[2])

                    # patient axis mappings: for DICOM->Nifti, depends purely on standards
                    out += '    "patient_axes_map": ["-x", "-y", "z"],\n'

                if self.use_exdcm:
                    if self.exdcm_path:
                        fn = e.file[sorted(e.file.keys())[0]]
                    else:
                        fn = basename(e.file[sorted(e.file.keys())[0]])
                    out += '    "exdcm": "%s",\n' % fn

                if filenames.has_key((study[0],study[1],k)):
                    fn = basename(filenames[study[0],study[1],k])
                    if fn.endswith('.nii') or fn.endswith('.nii.gz'):
                        out += '    "nii": "%s",\n' % fn
                    elif fn.endswith('.gipl') or fn.endswith('.gipl.gz'):
                        out += '    "gipl": "%s",\n' % fn

                out += '    "date": "%s",\n' % e.date
                out += '    "time": "%s"\n' % e.time
                out += '  }'
            out += '\n  ]\n'
            out += "}"
        out += ']\n'
        return(out)

if __name__ == "__main__":
    pass

# vim:et:sts=4:sw=4:tw=0
