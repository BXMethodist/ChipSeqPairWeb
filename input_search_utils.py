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
import re, pandas as pd
from collections import defaultdict
from difflib import SequenceMatcher


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

def spliterFinder(title, keyword):
    # find the spliter in the title, and return the keywords index and the spliter
    spliter = None
    index = None
    title = title.lower()

    space = title.count(" ")
    underscore = title.count("_")
    hyphen = title.count("-")
    choices = [" ", "_", "-"]
    counts = [space, underscore, hyphen]
    if space > 0 or underscore >0 or hyphen >0:
        max_count = 0
        choice = None
        for i in range(len(choices)):
            if counts[i] > max_count:
                max_count = counts[i]
                choice = i
        spliter = choices[choice]

    if spliter is None:
        if title.find(keyword.lower()) != -1:
            return None, -1
        else:
            return None, None
    elements = title.split(spliter)
    for i in range(len(elements)):
        if elements[i] == keyword.lower():
            index = i
            break
    return spliter, index


def Similarity(title1, keyword1, title2, keyword2, ignorecase):
    score = SequenceMatcher(None, title1, title2).ratio()

    if ignorecase:
        title1 = title1.lower().replace(keyword1.lower(), "").replace("chip-seq", "")
        title2 = title2.lower().replace(keyword2.lower(), "").replace("chip-seq", "")

        title1 = re.sub(r'rep[0-9]*', '', title1)
        title2 = re.sub(r'rep[0-9]*', '', title2)
    else:
        title1 = title1.replace(keyword1, "").lower().replace("chip-seq", "")
        title2 = title2.replace(keyword2, "").lower().replace("chip-seq", "")

        title1 = re.sub(r'rep[0-9]*', '', title1)
        title2 = re.sub(r'rep[0-9]*', '', title2)

    score_replace = SequenceMatcher(None, title1, title2).ratio()

    if score > score_replace:
        return score, score_replace
    else:
        return score_replace, score
    #return score_replace

def keyword(message, features, features_begin, ignorecase):
    if ignorecase:
        for feature in features:
            feature = feature.replace(" ","")
            feature = feature.replace("-", "")
            feature = feature.replace("_", "")
            if re.search(feature, message, flags=re.IGNORECASE):
                return feature
        for feature in features_begin:
            feature = feature.replace(" ","")
            feature = feature.replace("-", "")
            feature = feature.replace("_", "")
            if re.match(feature, message, flags=re.IGNORECASE):
                return feature
    else:
        for feature in features:
            feature = feature.replace(" ","")
            feature = feature.replace("-", "")
            feature = feature.replace("_", "")
            if re.search(feature, message):
                return feature
        for feature in features_begin:
            feature = feature.replace(" ","")
            feature = feature.replace("-", "")
            feature = feature.replace("_", "")
            if re.match(feature, message):
                return feature


def Character_Similarity(sample1, sample2):
    title1 = re.sub("\d", "", sample1.title)
    title2 = re.sub("\d", "", sample2.title)
    score = similar(title1, title2)
    for key, value in sample1.characteristics.items():
        if key in sample2.characteristics:
            score += similar(value, sample2.characteristics[key])

    # print sample1.id, sample2.id, score
    return score


def has_features(message, features, features_begin, ignorecase):
    if ignorecase:
        for feature in features:
            if re.search(feature, message, flags=re.IGNORECASE):
                return True
        for feature in features_begin:
            if re.match(feature, message, flags=re.IGNORECASE):
                return True
    else:
        for feature in features:
            if re.search(feature, message):
                return True
        for feature in features_begin:
            if re.match(feature, message):
                return True
    return False


def has_antibody(sample, keyword):
    for v in sample.antibody.values():
        if v.find(keyword) != -1:
            return True
    return False


def equal_antibody(sample, keyword):
    for v in sample.antibody.values():
        if v == keyword:
            return True
    return False


def isInput(sample, feature_key_word):
    non_capital_keywords = ['input','wce']

    if feature_key_word.find("H3K") == -1:
        capital_keywords = ['IgG']
    else:
        capital_keywords = ['IgG', '_H3_', " H3"]

    for c in non_capital_keywords:
        if sample.title.lower().find(c) != -1:
            return True, c
        if has_antibody(sample, c):
            return True, c

    for n in capital_keywords:
        if n ==' H3' and sample.title[-3:] == ' H3':
            return True, 'H3'
        elif sample.title.find(n) != -1 and n !=' H3':
            if n == "_H3_":
                return True, 'H3'
            return True, n
    return False, ""


def input_finder(output_surffix, output_path, HumanSamples, groupByGSE, encodeGSE, relatedSamples,
                 features, features_begin, ignorecase, output_type):
    FirstSampleToInput = defaultdict(set)

    ThirdSampleToInput = defaultdict(set)

    # get all the sample with key word in title
    titleCandidates = set()

    noneTitle = set()
    for key, value in HumanSamples.items():
        if value.title_found == True:
            titleCandidates.add(key)
        else:
            noneTitle.add(key)

    # print "title and none title", len(titleCandidates), len(noneTitle)
    not_found = 0
    # get their related samples
    for candidate in titleCandidates:
        sample = HumanSamples[candidate]

        feature_key_word = keyword(sample.title, features, features_begin, ignorecase)

        # print feature_key_word

        targetGSEs = set(sample.series)

        bestMatchID = set()
        bestSimilarity = float("-inf")
        bestSimilarity2 = float("-inf")
        input_keyword = ""

        for gse in targetGSEs:
            for relatedSample in groupByGSE[gse]:
                if relatedSample == candidate:
                    continue
                score = None
                boo, word = isInput(relatedSamples[relatedSample], feature_key_word)
                # print relatedSamples[relatedSample].title, boo, word
                if boo \
                        and sample.cellLine == relatedSamples[relatedSample].cellLine \
                        and sample.cellType == relatedSamples[relatedSample].cellType \
                        and sample.tissue == relatedSamples[relatedSample].tissue:

                    score1, score2 = Similarity(sample.title, feature_key_word, relatedSamples[relatedSample].title,
                                       word, ignorecase)

                    if score1 is not None and score1 > bestSimilarity:
                        bestSimilarity = score1
                        bestSimilarity2 = score2
                        bestMatchID = set()
                        bestMatchID.add(relatedSamples[relatedSample].id)
                        input_keyword = word
                    elif score1 is not None and score1 == bestSimilarity:
                        if input_keyword == 'H3' and word != 'H3':
                            bestMatchID = set()
                            bestMatchID.add(relatedSamples[relatedSample].id)
                            input_keyword = word
                            bestSimilarity = score1
                            bestSimilarity2 = score2
                        elif input_keyword != 'H3' and word == 'H3':
                            pass
                        else:
                            if score2 is not None and score2 > bestSimilarity2:
                                bestSimilarity2 = score2
                                bestMatchID = set()
                                bestMatchID.add(relatedSamples[relatedSample].id)

        if bestMatchID:
            FirstSampleToInput[sample.id] = FirstSampleToInput[sample.id].union(bestMatchID)
        else:
            not_found += 1

        # if candidate == 'GSM838673':
        #     print(FirstSampleToInput[sample.id])

    for key in noneTitle:
        sample = HumanSamples[key]
        targetGSEs = set(sample.series)

        best_char_score = 0
        best_id = set()

        for gse in targetGSEs:
            for relatedSample in groupByGSE[gse]:
                char_score = None
                for v in relatedSamples[relatedSample].antibody.values():
                    if v.find("input")!= -1 and sample.id != relatedSamples[relatedSample].id \
                            and sample.cellLine == relatedSamples[relatedSample].cellLine:
                        char_score = Character_Similarity(sample, relatedSamples[relatedSample])
                    elif v.lower().find(" wce ") != -1 and sample.id != relatedSamples[relatedSample].id \
                            and sample.cellLine == relatedSamples[relatedSample].cellLine:
                        char_score = Character_Similarity(sample, relatedSamples[relatedSample])
                    elif v.lower().find("whole cell extract") != -1 and sample.id != relatedSamples[relatedSample].id \
                            and sample.cellLine == relatedSamples[relatedSample].cellLine:
                        char_score = Character_Similarity(sample, relatedSamples[relatedSample])
                    elif v.find("IgG") != -1 and sample.id != relatedSamples[relatedSample].id \
                            and sample.cellLine == relatedSamples[relatedSample].cellLine:
                        char_score = Character_Similarity(sample, relatedSamples[relatedSample])
                if char_score is not None and char_score > best_char_score:
                    best_id = set()
                    best_id.add(relatedSamples[relatedSample].id)
                    best_char_score = char_score
                elif char_score is not None and char_score == best_char_score:
                    best_id.add(relatedSamples[relatedSample].id)

        if best_id:
            ThirdSampleToInput[sample.id] = best_id
        else:
            not_found+=1

    # print not_found
    output_type = output_type.replace(" ", "_")
    output1 = output_path + "First_" + output_surffix + "_" + output_type + "_Sample_To_Input.csv"
    output3 = output_path + "Third_" + output_surffix + "_" + output_type +"_Sample_To_Input.csv"

    table = []
    for key, value in FirstSampleToInput.items():
        row = [key] + [HumanSamples[key].title] + [HumanSamples[key].antibody]
        for id in value:
            row += [id] + [relatedSamples[id].title] + [relatedSamples[id].antibody]
        table.append(row)

    df = pd.DataFrame(table, columns=None)
    # df.to_csv(output1, sep=',', encoding='utf-8')

    table = []

    for key, value in ThirdSampleToInput.items():
        row = [key] + [HumanSamples[key].antibody]
        for id in value:
            row += [id] + [relatedSamples[id].antibody]
        table.append(row)
    df = pd.DataFrame(table, columns=None)
    # df.to_csv(output3, sep=',', encoding='utf-8')

    return FirstSampleToInput, ThirdSampleToInput



