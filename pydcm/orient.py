#!/usr/bin/env python
#
# Volconv - geometry-aware DICOM-to-NIfTI converter
# Oriented image class; representations are stored with DICOM conventions
# and functions are included to output in NIfTI conventions
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

from numpy import array, matrix, linalg
import numpy as np
import math
import re
# import pdb

# vector/product values below eps are considered to be zero
eps = 1e-5

# angular error allowed between normk and delta vectors (degrees)
atol = 2.0

def flip(s):
    """flip sign of a single-character axis label"""

    if re.match(r'-.', s):
        return s[1:]
    else:
        return '-' + s

def map_axis(s,axes):
    """
    Return an axis identified mapped via an axes permutation list:
    i,j,k are input axes; I,J,K are output axes.

    This is used to check how a given axis (eg phase direction) has
    been changed in the image grid by flips and tranpositions.  It
    operates on voxel indices, not world co-ordinates.
    """

    ret_axes = ['I','J','K']

    # exact match: "j" or "-j"
    if s in axes:
        return ret_axes[axes.index(s)]

    # "j" = s, "-j" in axes
    elif ("-" + s) in axes:
        return "-" + ret_axes[axes.index("-" + s)]

    # "-j" = s, "j" in axes
    elif re.match(r'-.', s) and s[1:] in axes:
        return "-" + ret_axes[axes.index(s[1:])]

    else:
        return "F"

class OrientedImage:

    def __init__(self,data,pixdim,orient,offset=[0.0,0.0,0.0],delta=None):
        """
        data = numpy array of image data on grid
        pixdim = mm dimensions of pixel sizes along each grid axis
        orient = [ix iy iz jx jy jz] orientation as DICOM [x y z] vectors
        offset = [x y z] corner voxel position of 0th slice
        delta = [x y z] vector between corner voxels of 0th and 1st slice

        DICOM defines the normk vector as (i x j): this vector triplet
        defines the image grid axes.  Unfortunately, depending on the
        sorting criteria used to stack the images, the actual delta
        (and thus positive indexing direction) may be along -normk.

        The way to ensure delta // normk is the --slice-3d argument.

        The convention in volconv is:

          * [i j k] always refers to [i j normk], ie. they refer to
            theoeretical grid-relative axes and are independent of the
            actual slice order within the volume.  

          * [x y z] always refer to DICOM world co-ordinates.
            
        This is the storage convention of all data structures in
        volconv, including this OrientedImage, and of all vector
        propreties output subsequently in index.json.  This makes
        index.json independent of the conversion process, and makes
        volconv somewhat agnostic about output formats.

          * If the stacking order is negative (delta = -normk), a
            flip is immediately recorded on k, so the image grid axes
            become [i j -k]

          * All subsequent updates to the image grid are applied as
            further transposition/negation to self.axes

        Mapping DICOM [x y z] to Nifti [X Y Z] co-ordinates is easy:
        [X Y Z] = [-x -y z].  The origin never changes.  This is a
        property of the specifications, not dependent on images.

        Mapping DICOM [i j normk] to Nifti [I J K] depends on the
        conversion which in turn depends on other fields.  But the
        self.axes field (which appears in index.json as "axes") will
        always hold enough information to do this conversion.
        """

        self.data = data
        self.pixdim = [float(x) for x in pixdim] # in mm

        if len(orient) == 1:
            self.i = [float(x) for x in orient[0][0:3]] # in DICOM cs
            self.j = [float(x) for x in orient[0][3:6]] # in DICOM cs
            self.mixed = False
        else:
            self.i = [1.0,0.0,0.0]
            self.j = [0.0,1.0,0.0]
            self.mixed = True

        self.offset = [float(x) for x in offset] # in DICOM cs
        self.delta = delta

        # original DICOM grid space
        self.axes = ['i','j','k']
        
        k = self.normk()
        qfac = self.checkSliceDir(k)

        if qfac < 0:
            self.axes[2] = '-k'

    def flipV(self):
        """v-flip; unit vectors remain in DICOM cs"""

        infovj = (self.data.shape[1]-1) * self.pixdim[1]
        for n in range(0,3):
            self.offset[n] = self.offset[n] + self.j[n] * infovj
        self.j = [-x for x in self.j]
        self.data = self.data[:,::-1,:]

        self.axes[1] = flip(self.axes[1])
    
    def flipH(self):
        """h-flip; unit vectors remain in DICOM cs"""

        infovi = (self.data.shape[0]-1) * self.pixdim[0]
        for n in range(0,3):
            self.offset[n] = self.offset[n] + self.i[n] * infovi
        self.i = [-x for x in self.i]
        self.data = self.data[::-1,:,:]

        self.axes[0] = flip(self.axes[0])
    
    def simplify(self,vector):
        """find a nearby 'simple' integer axis-aligned unit-vector"""

        largest = 0.0
        which   = -1

        for n in range(0,len(vector)):
            if abs(vector[n]) > abs(largest):
                largest = vector[n]
                which = n

        simplified = [0, 0, 0]
        if largest >= 0.0:
            simplified[which] = 1
        else:
            simplified[which] = -1
        
        return simplified

    def findOrient(self,style='long'):
        """find closest orthogonal orientation"""

        si = self.simplify(self.i)
        sj = self.simplify(self.j)

        if self.mixed:
            if style=='long': orient = "Mixed"
            else: orient = 'mix'

        elif si == [1,0,0] and sj == [0,1,0]:
            if style=='long': orient = "Axial"
            else: orient = 'axi'

        elif si == [0,1,0] and sj == [0,0,-1]:
            if style=='long': orient = "Sagittal"
            else: orient = 'sag'

        elif si == [1,0,0] and sj == [0,0,-1]:
            if style=='long': orient = "Coronal"
            else: orient = 'cor'

        else:
            if style=='long': orient = "Nonstd: %s, %s" % (repr(si), repr(sj))
            else: orient = 'obl'

        return orient

    def reOrient(self,new):
        """re-orient to given plane"""

        old  = self.findOrient()
        qfac = self.checkSliceDir()

        if old == new:
            return True

        if (old,new) == ("Coronal","Axial"):
            # transform: i'=i, j'=k, k'=-j
            self.data = self.data.swapaxes(1,2) # swap j and k
            self.j = self.normk()
            self.pixdim = [self.pixdim[0],self.pixdim[2],self.pixdim[1]]
            self.axes = [self.axes[0], self.axes[2], self.axes[1]]

            # flip k'
            self.data = self.data[:,:,::-1]
            infovk = (self.data.shape[2]-1) * self.pixdim[2]
            for n in range(0,3):
                self.offset[n] = self.offset[n] - self.normk()[n] * infovk
            self.axes[2] = flip(self.axes[2])

            # flip j' if k sign was wrong
            if qfac < 0:
                self.data = self.data[:,::-1,:] # j'=-k
                infovj = (self.data.shape[1]-1) * self.pixdim[1]
                for n in range(0,3):
                    self.offset[n] = self.offset[n] - self.j[n] * infovj
                self.axes[1] = flip(self.axes[1])

            self.recalcDelta()
            return True

        if (old,new) == ("Sagittal","Axial"):
            # transform: i'=-k, j'=i, k'=-j
            self.data = self.data.transpose(2,0,1)
            k = self.normk()
            self.j = self.i
            self.i = [-x for x in k]
            self.pixdim = [self.pixdim[2],self.pixdim[0],self.pixdim[1]]
            self.axes = [flip(self.axes[2]), self.axes[0], self.axes[1]]

            # flip k'
            self.data = self.data[:,:,::-1]
            infovk = (self.data.shape[2]-1) * self.pixdim[2]
            for n in range(0,3):
                self.offset[n] = self.offset[n] - self.normk()[n] * infovk
            self.axes[2] = flip(self.axes[2])
            
            # flip i' (again) if k sign was correct
            if qfac > 0:
                self.data = self.data[::-1,:,:] # i'=-l
                infovi = (self.data.shape[0]-1) * self.pixdim[0]
                for n in range(0,3):
                    self.offset[n] = self.offset[n] - self.i[n] * infovi
                self.axes[0] = flip(self.axes[0])

            self.recalcDelta()
            return True
        
        return False

    def qdata(self):
        """returns qdata in NIFTI cs"""
        q = self.quaternion()
        qfac = q[0]
        q1 = q[2:5]
        q1.extend([-self.offset[0], -self.offset[1], self.offset[2]])
        return (qfac, q1)

    def normk(self):
        """return k as positive cross-product i x j"""
        k = [0.0, 0.0, 0.0]
        i = self.i
        j = self.j
        k[0] = i[1]*j[2] - i[2]*j[1]
        k[1] = i[2]*j[0] - i[0]*j[2]
        k[2] = i[0]*j[1] - i[1]*j[0]
        return k

    def recalcDelta(self):
        """recalculate delta, assuming grid is orthogonal"""
        k = self.normk()
        s = self.pixdim[2]
        self.delta = [ki * s for ki in k]
   
    def useSliceGap(self):
        """calculate slice gap from delta"""
        answer = 0

        if self.delta == None:
            return 1.0

        magn = math.sqrt((self.delta[0] ** 2.) +
            (self.delta[1] ** 2.) +
            (self.delta[2] ** 2.))
        self.pixdim[2] = magn
        return magn

    def checkSliceDir(self,k=None):
        """
        checks whether k really points along slice-dir, return sign.

         1 means normk ==   delta / |delta|
        -1 means normk == - delta / |delta|
         0 means delta is shorter than eps
        """

        # print "\nNormk:  ", self.normk()
        # print "Normk:  ", self.normk()
        # print "Pixdim: ", self.pixdim
        # print "Delta:  ", self.delta

        if k == None:
            k = self.normk()

        if self.delta == None:
            return 1.0

        if linalg.norm(self.delta) < eps:
            return 0.0

        kv = np.array(k)
        dv = np.array(self.delta)
        normdot = np.dot(kv,dv)/(linalg.norm(kv)*linalg.norm(dv))

        if normdot > 1.0 and normdot < 1.0+eps:
            normdot = 1.0
        
        if normdot < -1.0 and normdot > -1.0-eps:
            normdot = -1.0

        angle = np.arccos(normdot)*180.0/np.pi
        
        # print "kv", kv
        # print "dv", dv
        # print "kv norm", linalg.norm(kv)
        # print "dv norm", linalg.norm(dv)
        # print "dot", np.dot(kv,dv)
        # print "norm dot", normdot
        # print "angle", angle

        if angle > (180.0-atol):
            return -1.0
        elif angle < atol:
            return 1.0
        else:
            raise ValueError, ("interslice vector and normal vector differ " + \
                "by > %.2f deg (actual angle: %.2f deg)" % (atol,angle,))

    def quaternion(self):
        """calculate quaternion repr from NIFTI C library; in NIFTI cs"""

        r11 = - self.i[0] # negation DICOM->NIFTI
        r21 = - self.i[1] # negation DICOM->NIFTI
        r31 = self.i[2]
        
        r12 = - self.j[0] # negation DICOM->NIFTI
        r22 = - self.j[1] # negation DICOM->NIFTI
        r32 = self.j[2]
      
        # compute k in DICOM cs, check sign of qfac
        k = self.normk()
        qfac = self.checkSliceDir(k)

        # now find columns in NIFTI cs (negate again)
        r13 = - k[0] # negation DICOM->NIFTI
        r23 = - k[1] # negation DICOM->NIFTI
        r33 = k[2]

        a = r11 + r22 + r33 + 1.0

        if a > 0.5:
            a = 0.5 * math.sqrt(a)
            b = 0.25 * (r32-r23) / a
            c = 0.25 * (r13-r31) / a
            d = 0.25 * (r21-r12) / a
        else:
            xd = 1.0 + r11 - (r22+r33)
            yd = 1.0 + r22 - (r11+r33)
            zd = 1.0 + r33 - (r11+r22)
            if xd > 1.0:
                b = 0.5 * math.sqrt(xd)
                c = 0.25* (r12+r21) / b
                d = 0.25* (r13+r31) / b
                a = 0.25* (r32-r23) / b
            elif yd > 1.0:
                c = 0.5 * math.sqrt(yd)
                b = 0.25* (r12+r21) / c
                d = 0.25* (r23+r32) / c
                a = 0.25* (r13-r31) / c
            else:
                d = 0.5 * math.sqrt(zd)
                b = 0.25* (r13+r31) / d
                c = 0.25* (r23+r32) / d
                a = 0.25* (r21-r12) / d
        if a < 0.0:
            b=-b
            c=-c
            d=-d
            a=-a
        
        return [qfac,a,b,c,d]

    def dcm_to_grid(self, v):
        """convert a vector from DICOM co-ordinates to the DICOM
        grid co-ordinates [i j normk]; note that we do not care about
        the actual stacking direction here."""

        # only use on un-transposed image orientation
        # - don't care about k, because we're using self.normk()
        assert self.axes[0] == "i"
        assert self.axes[1] == "j"

        k = self.normk()

        T = matrix([self.i, self.j, k]).transpose()
        Ti = linalg.inv(T)

        vg = Ti*matrix(v).transpose()
        vgl = array(vg)[:,0].tolist()
        return vgl

    def map_axis(self,s):
        """apply map_axis to show where a particular original DICOM
        axes (i,j,k) has ended up in the permuted version of this
        image"""

        return map_axis(s,self.axes)
