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

import re, sqlite3, pandas as pd, os, json
from collections import defaultdict
from sample import GSM
from input_search_utils import input_finder
from encode import encode_search
import re

def Search(output_prefix, output_path, keywords, _Species, cwd,
           _CellLines, _CellTypes, _Organs, _Tissues, _Sequencing_protocol='chip-seq', _inputEmail=None):

    chip_db = sqlite3.connect(cwd)
    df = pd.read_sql_query('SELECT * from metadata', con=chip_db, index_col=['Data_ID'])

    df = df[df['Species'].str.lower() == _Species.lower()]
    GSEGSM_map = defaultdict(set)

    antibody_df = df.copy()

    avals_lowers = [x.lower() for x in os.listdir('./aval')]
    avals = ['./aval/' + x for x in os.listdir('./aval')]
    keywords_with_species = [word.lower()+'_'+_Species.replace(' ', '').lower()+'.csv' for word in keywords]

    avals_results = []

    for i in range(len(avals_lowers)):
        if avals_lowers[i] in keywords_with_species:
            avals_results.append(i)

    if len(avals_results) != 0:
        curated_result_df = pd.read_csv(avals[avals_results[0]], index_col=0)
        curated_index = curated_result_df.index
    else:
        curated_index = None

    keywords = [word.lower() for word in keywords]

    pattern = '|'.join(map(re.escape, keywords))

    inputs_pattern = '|'.join(map(re.escape, ['input', 'control', 'igg', 'wce']))

    df = df[(df.Title.str.lower().str.contains(pattern)) &
                       (~((df.Title.str.lower().str.contains(inputs_pattern)) |
                          (df.Antibody.str.lower().str.contains(inputs_pattern))))]

    antibody_df = df[(df.Antibody.str.lower().str.contains(pattern)) &
                       (~((df.Title.str.lower().str.contains(inputs_pattern)) |
                          (df.Antibody.str.lower().str.contains(inputs_pattern))))]


    celllines = [word.lower() for word in _CellLines]
    cellline_pattern = '|'.join(map(re.escape, celllines))
    df = df[df.CellLine.str.lower().str.contains(cellline_pattern)]

    antibody_df = antibody_df[antibody_df.CellLine.str.lower().str.contains(cellline_pattern)]

    celltypes = [word.lower() for word in _CellTypes]
    celltype_pattern = '|'.join(map(re.escape, celltypes))
    df = df[df.CellType.str.lower().str.contains(celltype_pattern)]
    antibody_df = antibody_df[antibody_df.CellType.str.lower().str.contains(celltype_pattern)]

    organs = [word.lower() for word in _Organs]
    organ_pattern = '|'.join(map(re.escape, organs))
    df = df[df.Organ.str.lower().str.contains(organ_pattern)]
    antibody_df = antibody_df[antibody_df.Organ.str.lower().str.contains(organ_pattern)]

    tissues = [word.lower() for word in _Tissues]
    tissue_pattern = '|'.join(map(re.escape, tissues))
    df = df[df.Tissue.str.lower().str.contains(tissue_pattern)]
    antibody_df = antibody_df[antibody_df.Tissue.str.lower().str.contains(tissue_pattern)]

    protocols = [word.lower() for word in _Sequencing_protocol]
    protocol_pattern = '|'.join(map(re.escape, protocols))
    df = df[df.SequencingProtocol.str.lower().str.contains(protocol_pattern)]
    antibody_df = antibody_df[antibody_df.SequencingProtocol.str.lower().str.contains(protocol_pattern)]

    title_indexes = df.index
    antibody_indexes = antibody_df.index

    confidences = []

    Data_desciption = []

    for id in df.index:
        gses = df.ix[id, 'Study_ID'].split(',')
        for gse in gses:
            GSEGSM_map[gse].add(id)

        if id in title_indexes and id in antibody_indexes:
            confidences.append('High Confidence')
            desciption = 'Title:'+df.ix[id, 'Title'] +';'+'Antibody:'+ df.ix[id, 'Antibody']
            Data_desciption.append(desciption)
        elif id in title_indexes:
            confidences.append('Medium Confidence')
            desciption = 'Title:' + df.ix[id, 'Title']
            Data_desciption.append(desciption)
        elif id in antibody_indexes:
            confidences.append('Low Confidence')
            desciption = 'Antibody:' + df.ix[id, 'Antibody']
            Data_desciption.append(desciption)
        else:
            confidences.append('No Confidence')
            desciption = ''
            Data_desciption.append(desciption)

    df['Confidence'] = confidences
    df['Data Description'] = Data_desciption

    confidence_pattern = '|'.join(map(re.escape, ['High', 'Medium', 'Low']))

    result_df = df[df.Confidence.str.contains(confidence_pattern)].copy()

    samples = dfToSamples(result_df)

    inputs = []
    inputs_descriptions = []

    if keywords != []:
        Category1, Category3, relatedSamples = Input(output_prefix, output_path, keywords, _Species, GSEGSM_map, samples, df)
        for sample_id in result_df.index:
            input_id = ''
            cur_description = ''
            if sample_id in Category1:
                input_ids = Category1[sample_id]
                input_id = input_id + ','.join(list(input_ids))

                for id in input_ids:
                    cur_title = relatedSamples[id].title
                    cur_description += 'Title:' + cur_title + '; '
                    cur_antibodies = json.dumps(relatedSamples[id].antibody)
                    cur_description += 'Antibody:'+ cur_antibodies

            elif sample_id in Category3:
                input_ids = Category3[sample_id]
                input_id = input_id + ','.join(list(input_ids))

                for id in input_ids:
                    cur_title = relatedSamples[id].title
                    cur_description += 'Title:' + cur_title + '; '
                    cur_antibodies = json.dumps(relatedSamples[id].antibody)
                    cur_description += 'Antibody:' + cur_antibodies
            inputs.append(input_id)
            inputs_descriptions.append(cur_description)
    else:
        for sample_id in result_df.index:
            input_id = ''
            cur_description = ''
            inputs.append(input_id)
            inputs_descriptions.append(cur_description)


    result_df['Input'] = inputs
    result_df['Input Description'] = inputs_descriptions

    result_df = result_df.ix[:, ['Study_ID', 'Confidence','Data Description', 'Title', 'Antibody', 'Input',
                                 'Input Description', 'InstrumentModel', 'RawData', 'SequencingProtocol',
                                 'Species', 'CellLine', 'CellType',
                                  'Tissue', 'Organ']]

    result_df.columns = ["Study_ID", 'Confidence', "Data_Description", "Title", "Experiment target/antibody",
               "Input", "Input_Description", "Instrument_Model", "Raw Data", "Sequencing_Protocol",
               "Species", "Cell Line", "Cell Type", 'Tissue', 'Organ']

    samples_encode, human_encode, human_encode_map = encode_search(output_prefix, keywords,
                                                                   _CellLines, _CellTypes,
                                                                   _Organs, _Tissues,
                                                                   _Sequencing_protocol, output_type=_Species)

    if human_encode is not None:
        result_df = result_df.append(human_encode)
        samples.update(human_encode_map)

    # output_name = output_path + output_prefix + '_' + '_'.join(_Species.split()) + '.csv'
    # result_df.to_csv(output_name)
    if curated_index is not None:
        result_df = result_df.ix[curated_index, :]
        for sample_id in samples.keys():
            if sample_id not in curated_index:
                del samples[sample_id]

    return result_df, samples

def dfToSamples(df):
    samples = {}
    for id in df.index:
        samples[id] = PandasToSample(id, df.ix[id, :])
    return samples

def PandasToSample(id, sample_series):
    sample=GSM(id)
    sample.series = sample_series['Study_ID'].split(',')
    sample.characteristics = json.loads(sample_series['Characteristics'])
    sample.libraryStrategy = sample_series['SequencingProtocol']
    sample.SRA = sample_series['RawData']
    sample.features = sample_series['Data Description']
    sample.title = sample_series['Title']
    sample.InstrumentID = sample_series['InstrumentModel']
    sample.organism = sample_series['Species']
    sample.antibody = json.loads(sample_series['Antibody'])
    sample.tissue = sample_series['Tissue']
    sample.cellLine = sample_series['CellLine']
    sample.cellType = sample_series['CellType']

    if sample_series['Confidence'] == 'High Confidence':
        sample.title_ab = True
        sample.title_found = True
        sample.ab_found = True
    elif sample_series['Confidence'] == 'Medium Confidence':
        sample.title_found = True
    elif sample_series['Confidence'] == 'Low Confidence':
        sample.ab_found = True

    sample.organ = sample_series['Organ']
    return sample


def Input(output_prefix, output_path,features, speices, GSEGSM_map, samples, df):
    encodeGSE = set()
    features_begin = []
    ignorecase = True

    relatedSamples = {}

    for sample_id in samples.keys():
        cur_relatedGSEs = samples[sample_id].series

        for gse in cur_relatedGSEs:
            cur_relatedGSMs = GSEGSM_map[gse]

            for gsm in cur_relatedGSMs:
                relatedSamples[gsm] = PandasToSample(gsm, df.ix[gsm, :])

    Category1, Category3 = input_finder(output_prefix, output_path, samples, GSEGSM_map, encodeGSE, relatedSamples,
                                        features, features_begin, ignorecase, speices)
    return Category1, Category3, relatedSamples


if __name__ == "__main__":
    results = Search('' ,'',[],'Homo sapiens','/Users/boxia/PycharmProjects/ChiPSeqPairWeb/db/chipseq.db',['mcf-7'],[],[],[])
    print results
