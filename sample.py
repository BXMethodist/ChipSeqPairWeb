"""
Copyright (c) <2017> <Dr. Kaifu Chen lab, Research Institute, Houston Methodist Hospital >

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import re
from urllib.request import urlopen

def search_term_to_GSM(terms):
    result_ids = set()
    result_gsms = set()

    for term in terms:
        cur_result_ids = readid("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSM%5BETYP%5D"+term+"&retmax=10000")

        result_ids = result_ids.union(cur_result_ids)

    result_ids = list(result_ids)

    for i in range(0, len(result_ids), 500):
        ids = ",".join(result_ids[i:i+500])
        gsm_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&version=2.0&id=" + ids

        f = urlopen(gsm_url)
        content = str(f.read(), 'utf-8')
        cur_gsms = [x.strip() for x in re.findall('\sGSM[0-9]+\s', content)]
        result_gsms = result_gsms.union(cur_gsms)
    gsms = set()

    for gsm in result_gsms:
        gsms.add(gsm.strip())

    return gsms


def readid(url):
    f = urlopen(url)
    content = str(f.read())
    content = content[content.find("<IdList>") + 8:content.find("</IdList>")-1]
    content = content.replace("<Id>", " ")
    content = content.replace("</Id>", " ")
    result_ids = set(content.split())
    return result_ids


class GSM:
    def __init__(self, GSM_id):
        self.id = GSM_id
        self.series = []
        self.platForm = ""
        self.characteristics = {}
        self.supplementaryData = {}
        self.relations = {}
        self.libraryStrategy = ""
        self.SRA = ""
        self.type = ""
        self.features = ""
        self.title = ""
        self.InstrumentID = ""
        self.organism = ""
        self.url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="+GSM_id+"&targ=self&form=text&view=quick"

    # potential attributes
        self.antibody = {}
        self.treatment = {}
        self.tissue = ""
        self.disease = ""
        self.cellLine = ""
        self.cellType = ""
        self.genotype = {}
        self.title_found = False
        self.ab_found = False
        self.title_ab = False
        self.input = ""
        self.encode = ""
        self.roadmap = ""
        self.organ = ""
