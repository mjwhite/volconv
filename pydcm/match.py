#!/usr/bin/env python
#
# Volconv - geometry-aware DICOM-to-NIfTI converter
# Series matching parser
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

import ConfigParser
from dicom import fixser
from datetime import date
import dicom
import re

def dateDiff(date1,date2):
    """
    find the difference between two dates
    """
    d1 = date(int(date1[0:4]), int(date1[4:6]), int(date1[6:8]))
    d2 = date(int(date2[0:4]), int(date2[4:6]), int(date2[6:8]))
    delta = d2 - d1
    return delta.days

def parseRange(text):
    """
    convert a range string like 2-3, 10-20, 4-, -9, or 2 to a list
    containing the endpoints.  A missing endpoint is set to None.
    """
    def toNumeric(elt):
        if elt == "":
            return None
        else:
            return int(elt)

    if re.search(r'-', text):
        rng = text.split(r':')
        rng = [toNumeric(elt) for elt in rng]
    else:
        v = int(text)
        rng = [v, v]

    return rng

class NameMatcher:

    def __init__(self,file):
        self.conf = ConfigParser.RawConfigParser()
        self.conf.read(file)
        self.default_template =  "%(alias)?(-count)?(-t)?(-echo)"
        self.default_ignorecase = 1
        self.default_tidy = 1

        self._buildPatterns()

    def _buildPatterns(self):
        """
        build pattern dictionary from configuration file
        """

        self.patterns = {}
        aliases = self.conf.sections()

        if self.conf.has_section("default"):
            aliases.remove("default")

            if self.conf.has_option("default","template"):
                self.default_template = self.conf.get("default","template")

            if self.conf.has_option("default","ignorecase"):
                self.default_ignorecase = self.conf.getint("default","ignorecase")

            if self.conf.has_option("default","tidy"):
                self.default_tidy = self.conf.getint("default","tidy")

        for a in aliases:
            dic = {}
            self.patterns[a] = dic

            if self.conf.has_option(a,"pattern"):
                dic["pattern"] = self.conf.get(a,"pattern")

            if self.conf.has_option(a,"days"):
                dic["days"] = parseRange(self.conf.get(a,"days"))
            
            if self.conf.has_option(a,"type"):
                dic["type"] = self.conf.get(a,"type")

            if self.conf.has_option(a,"count"):
                dic["count"] = parseRange(self.conf.get(a,"count"))

            if self.conf.has_option(a,"series"):
                dic["series"] = parseRange(self.conf.get(a,"series"))
            
            if self.conf.has_option(a,"study"):
                dic["study"] = parseRange(self.conf.get(a,"study"))
            
            if self.conf.has_option(a,"template"):
                dic["template"] = self.conf.get(a,"template")
            
            if self.conf.has_option(a,"ignorecase"):
                dic["ignorecase"] = self.conf.getint(a,"ignorecase")
            
            if self.conf.has_option(a,"tidy"):
                dic["tidy"] = self.conf.getint(a,"tidy")

            dic["matched"] = 0   # how many have matched?
            dic["counted"] = 0   # how many matches were counted?

    def ignorecaseFlag(self,a):
        """
        get the ignorecase flag (or the default) for an alias
        """

        dic = self.patterns[a]
        if dic.has_key("ignorecase"):
            ignorecase = dic["ignorecase"]
        else:
            ignorecase = self.default_ignorecase

        if ignorecase:
            return re.IGNORECASE
        else:
            return 0
       
    def tidyIfRequired(self,a,name):
        """
        tidy a protocol name if flagged to do so for this match
        """

        dic = self.patterns[a]
        if dic.has_key("tidy"):
            tidy = dic["tidy"]
        else:
            tidy = self.default_tidy

        if tidy:
            return dicom.tidy_protoname(name)
        else:
            return name

    def findMatches(self,dcm):
        """
        find and store all matches in the given set of DICOM data (use
        the match() method to retrieve individual ones later)
        """

        series_list = []

        for study_id in dcm.studies.keys():
            study = dcm.studies[study_id]
            for series_id in dcm.studies[study_id].keys():
                series = study[series_id]
                desc = series.desc

                identifier = {
                        "study_no":   study_id[0],
                        "study_name": study_id[1],
                        "series_no":  series_id,
                        "series_int": int(re.sub(r'[^0-9]','',series_id)),
                        "date": series.date,
                        "time": series.time,
                        "stdate": series.stdate,
                        "sttime": series.sttime,
                        "desc": series.desc,
                        "type": series.imtype,
                }

                series_list.append(identifier)

        def order(a,b):
            ast = a["stdate"] + a["sttime"]
            bst = b["stdate"] + b["sttime"]
            date_cmp = cmp(ast,bst)
            if date_cmp == 0:
                return cmp(fixser(a["series_no"]), fixser(b["series_no"]))
            else:
                return date_cmp

        series_list.sort(order)

        baseline_date = series_list[0]["stdate"]
        self.matches = {}

        last_study = (series_list[0]["study_no"], series_list[0]["study_name"])
        study_count = 0

        for e in series_list:

            # increment the study count if necessary
            this_study = (e["study_no"], e["study_name"])
            if  this_study != last_study:
                study_count += 1
                last_study   = this_study

            for a in self.patterns.keys():
                p = self.patterns[a]

                if p.has_key("pattern"):
                    if not re.search(p["pattern"],
                            self.tidyIfRequired(a,e["desc"]),
                            self.ignorecaseFlag(a)):
                        continue
                
                if p.has_key("type"):
                    if not re.search(p["type"], e["type"],
                            self.ignorecaseFlag(a)):
                        continue

                if p.has_key("days"):
                    rng = p["days"]
                    age = dateDiff(baseline_date,e["stdate"])
                    if not (rng[0] is None) and age < rng[0]:
                        continue
                    if not (rng[1] is None) and age > rng[0]:
                        continue
                
                if p.has_key("study"):
                    rng = p["study"]
                    if not (rng[0] is None) and study_count < rng[0]:
                        continue
                    if not (rng[1] is None) and study_count > rng[0]:
                        continue
                
                # this is a match, but is it in the count range?
                true_count = p["matched"]
                offset_count = true_count
                p["matched"] += 1

                if p.has_key("count"):
                    rng = p["count"]
                    if not (rng[0] is None) and true_count < rng[0]:
                        continue
                    if not (rng[1] is None) and true_count > rng[1]:
                        continue
                    if not (rng[0] is None):
                        offset_count -= rng[0]

                if p.has_key("series"):
                    rng = p["series"]
                    series = e["series_int"]
                    if not (rng[0] is None) and series < rng[0]:
                        continue
                    if not (rng[1] is None) and series > rng[0]:
                        continue

                # save the match
                self.matches[(
                    e["study_no"],
                    e["study_name"],
                    e["series_no"],
                    )] = (a, offset_count)

                p["counted"] += 1

                # each series can only match one pattern
                break
    
    def template(self,study_no,study_name,series_no):
        """
        retrieve the template for this alias (or the default one)
        """

        key = (study_no, study_name, series_no)
        tpl = self.default_template

        if self.matches.has_key(key):
            alias, count = self.matches[key]
            p = self.patterns[alias]
            if p.has_key("template"):
                tpl = p["template"]

        return tpl

    def match(self,study_no,study_name,series_no):
        """
        retrieve the alias for a given series -- returns:
            None if there's no match
            (alias, None) if there's a unique match
            (alias, count) if there's a non-unique match
        """

        key = (study_no, study_name, series_no)
        if self.matches.has_key(key):
            alias, count = self.matches[key]
            if self.patterns[alias]["counted"] <= 1:
                count = None
            return (alias, count)
        else:
            return None
