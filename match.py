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

from collections import defaultdict
from input_search_utils import keyword, Similarity, Character_Similarity
import pandas as pd
from search import Search


def Match(output_prefix1, output_prefix2, output_path, first_features, second_features, species, cwd,
          _CellLines, _CellTypes, _Organs, _Tissues, _inputEmail):

    first_df, first_samples = Search(output_prefix1, output_path, first_features, species, cwd,
                                     _CellLines, _CellTypes, _Organs, _Tissues, _inputEmail=_inputEmail)
    second_df, second_samples = Search(output_prefix2, output_path, second_features, species, cwd,
                                       _CellLines, _CellTypes, _Organs, _Tissues, _inputEmail=_inputEmail)

    first_features_begin = []
    first_ignorecase = True
    second_begin_features = []
    second_ignorecase = True

    ## STEP1 groupbygse
    first_groupbygse = defaultdict(set)

    for sample in first_samples.values():
        for gse in sample.series:
            first_groupbygse[gse].add(sample.id)

    second_groupbygse = defaultdict(set)

    for sample in second_samples.values():
        for gse in sample.series:
            second_groupbygse[gse].add(sample.id)

    first_pairs = defaultdict(set)
    second_pairs = defaultdict(set)

    for sample in first_samples.values():
        targetGSMs = set()
        for gse in sample.series:
            targetGSMs = targetGSMs.union(second_groupbygse[gse])
        if sample.title_found:
            feature_key_word = keyword(sample.title, first_features, first_features_begin, first_ignorecase)

            best_score = float("-inf")
            best_id = set()

            for related_id in targetGSMs:
                score = None
                related_sample = second_samples[related_id]
                if related_sample.title_found:
                    related_feature_key_word = keyword(related_sample.title, second_features,
                                                       second_begin_features, second_ignorecase)

                    if sample.cellType == related_sample.cellType:
                        boo = (first_ignorecase or second_ignorecase)
                        score = Similarity(sample.title, feature_key_word, related_sample.title, related_feature_key_word,
                                           boo)[0]

                    if score is not None and score > best_score:
                        best_score = score
                        best_id = set()
                        best_id.add(related_sample.id)
                    elif score is not None and score == best_score:
                        best_id.add(related_sample.id)
                else:
                    continue

            if best_id:
                first_pairs[sample.id] = first_pairs[sample.id].union(best_id)

        else:
            # TODO none-title sample
            best_score = float("-inf")
            best_id = set()
            for related_id in targetGSMs:
                score = None
                related_sample = second_samples[related_id]

                if sample.cellType == related_sample.cellType:
                    score = Character_Similarity(sample, related_sample)

                if score is not None and score > best_score:
                    best_score = score
                    best_id = set()
                    best_id.add(related_sample.id)
                elif score is not None and score == best_score:
                    best_id.add(related_sample.id)
            if best_id:
                first_pairs[sample.id] = first_pairs[sample.id].union(best_id)

    for sample in second_samples.values():
        targetGSMs = set()
        for gse in sample.series:
            targetGSMs = targetGSMs.union(first_groupbygse[gse])
        if sample.title_found:
            feature_key_word = keyword(sample.title, second_features, second_begin_features, second_ignorecase)

            best_score = float("-inf")
            best_id = set()

            for related_id in targetGSMs:
                score = None
                related_sample = first_samples[related_id]
                if related_sample.title_found:
                    related_feature_key_word = keyword(related_sample.title, first_features,
                                                       first_features_begin, first_ignorecase)

                    if sample.cellType == related_sample.cellType:
                        boo = (first_ignorecase or second_ignorecase)
                        score = Similarity(sample.title, feature_key_word, related_sample.title, related_feature_key_word,
                                           boo)[0]

                    if score is not None and score > best_score:
                        best_score = score
                        best_id = set()
                        best_id.add(related_sample.id)
                    elif score is not None and score == best_score:
                        best_id.add(related_sample.id)
                else:
                    continue

            if best_id:
                second_pairs[sample.id] = second_pairs[sample.id].union(best_id)

        else:
            # TODO none-title sample
            best_score = float("-inf")
            best_id = set()
            for related_id in targetGSMs:
                score = None
                related_sample = first_samples[related_id]

                if sample.cellType == related_sample.cellType:
                    score = Character_Similarity(sample, related_sample)

                if score is not None and score > best_score:
                    best_score = score
                    best_id = set()
                    best_id.add(related_sample.id)
                elif score is not None and score == best_score:
                    best_id.add(related_sample.id)
            if best_id:
                second_pairs[sample.id] = second_pairs[sample.id].union(best_id)


    wrong_pairs = set()
    new_pairs = set()
    ### cross validation
    for key, values in first_pairs.items():
        for value in values:
            if key not in second_pairs[value]:
                wrong_pairs.add((key, value))
            elif (key, value) not in new_pairs:
                new_pairs.add((key, value))

    #### save file
    if not output_path.endswith("/"):
        output_path += "/"

    output = open(output_path + output_prefix1 + "_vs_" + output_prefix2 + ".csv", "w")

    table = []
    headers = ["Data1_ID", "Data1_Title", "Data1_descriptions", "Data1_Study_ID",
               "Data1_CellLine", "Data1_CellType", "Data1_Tissue",
               "Data2_id", "Data2_Title", "Data2_descriptions", "Data2_Study_ID",
               "Data2_CellLine", "Data2_CellType", "Data2_Tissue"]

    for pair in new_pairs:
        key, value = pair
        table.append([key, first_samples[key].title,
                      first_samples[key].features, first_samples[key].series,
                      first_samples[key].cellLine, first_samples[key].cellType, first_samples[key].tissue,
                      value, second_samples[value].title,
                      second_samples[value].features, second_samples[value].series,
                      second_samples[value].cellLine, second_samples[value].cellType, second_samples[value].tissue])

    df = pd.DataFrame(table, columns=headers)
    # df.to_csv(output, sep=',', encoding='utf-8', index=False)

    return first_df, second_df, df
