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

import pandas as pd
from sample import GSM

def encode_search(output_prefix, keywords, keywords_begin=(), type_seq='chip-seq',
                  candidates=(), ignorecase=True, output_type='Homo sapiens'):

    seq_types = ['RNA-PET', 'ATAC-seq', 'genotyping by high throughput sequencing assay', 'HiC', 'RNA-seq',
                 'genetic modification followed by DNase-seq', 'ChIP-seq', '5C', 'siRNA knockdown followed by RNA-seq',
                 'eCLIP', 'TAB-seq', 'CRISPR genome editing followed by RNA-seq', 'RRBS', 'RIP-seq',
                 'shRNA knockdown followed by RNA-seq', 'iCLIP', 'RNA Bind-n-Seq',
                 'single cell isolation followed by RNA-seq', 'RAMPAGE', 'FAIRE-seq', 'DNase-seq',
                 'CAGE', 'PAS-seq', 'whole-genome shotgun bisulfite sequencing', "5' RLM RACE", 'Repli-seq',
                 'ChIA-PET', 'DNA-PET', 'microRNA-seq']
    seq_types_map = {}
    for type in seq_types:
        seq_types_map[type.lower()] = type

    if type_seq.lower() in seq_types_map:
        url = 'https://www.encodeproject.org/metadata/type=Experiment&files.file_type=fastq&assay_title='+\
              seq_types_map[type_seq.lower()]+'/metadata.tsv'
    else:
        url = 'https://www.encodeproject.org/metadata/type=Experiment&files.file_type=fastq/metadata.tsv'

    human_encode_map = {}

    df = pd.read_csv(url, sep='\t', dtype=str)
    case = False if ignorecase else True

    features = "|".join(keywords)
    df = df[df['Experiment target'].str.contains(features, case=case, na=False)]
    if ignorecase:
        features_begin = "|".join(keywords_begin).lower()
        df = df[df['Experiment target'].str.lower().str.startswith(features_begin, na=False)]
    else:
        features_begin = "|".join(keywords_begin)
        df = df[df['Experiment target'].str.startswith(features_begin, na=False)]

    df['Confidence'] = ['High Confident'] * len(df.index)
    df['Input_Description'] = ['indicated by encode'] * len(df.index)
    df[output_prefix.capitalize() + "_Description"] = df['Experiment target']

    samples_df = df.copy()
    samples_df = samples_df[['File accession', 'Biosample term id', 'Assay', 'Biosample term name',
                             'Biosample type', 'Biosample organism', 'File download URL', 'Platform',
                             'Experiment target', output_prefix.capitalize() + "_Description",
                             'Confidence']]
    samples_df.columns = ['Data_ID', 'Study_ID', 'Sequencing_Protocol', 'Cell Line', 'Cell Type', 'Species',
                          'Raw Data', 'Instrument_Model', 'Experiment target/antibody',
                          output_prefix.capitalize() + "_Description", 'Confidence']
    samples_df = samples_df.set_index(['Data_ID'])

    df = df[df['Biosample organism'].str.contains(output_type, case=case, na=False)]
    df = df[df['Assay'].str.contains(type_seq, case=case, na=False)]

    if len(candidates) > 0:
        df = df.ix[candidates, ]
    elif candidates is None:
        return None, None, None

    for i in df.index:
        sample = GSM(df.ix[i, 'File accession'])
        sample.series = df.ix[i, 'Biosample term id'] if not pd.isnull(df.ix[i, 'Biosample term id']) else ""
        sample.libraryStrategy = df.ix[i, 'Assay'] if not pd.isnull(df.ix[i, 'Assay']) else ""
        sample.cellLine = df.ix[i, 'Biosample term name'] if not pd.isnull(df.ix[i, 'Biosample term name']) else ""
        sample.cellType = df.ix[i, 'Biosample type'] if not pd.isnull(df.ix[i, 'Biosample type']) else ""
        sample.organism = df.ix[i, 'Biosample organism'] if not pd.isnull(df.ix[i, 'Biosample organism']) else ""
        sample.antibody = df.ix[i, 'Experiment target']
        sample.features = df.ix[i, 'Experiment target']
        char = {}
        char['Biosample term name'] = df.ix[i, 'Biosample term name'] \
            if not pd.isnull(df.ix[i, 'Biosample term name']) else ""
        char['Biosample type'] = df.ix[i, 'Biosample type'] \
            if not pd.isnull(df.ix[i, 'Biosample type']) else ""
        char['Biosample life stage'] = df.ix[i, 'Biosample life stage'] \
            if not pd.isnull(df.ix[i, 'Biosample life stage']) else ""
        char['Biosample treatments'] = df.ix[i, 'Biosample treatments'] \
            if not pd.isnull(df.ix[i, 'Biosample treatments']) else ""
        char['Biosample subcellular fraction term name'] = df.ix[i, 'Biosample subcellular fraction term name'] \
            if not pd.isnull(df.ix[i, 'Biosample subcellular fraction term name']) else ""
        char['Biosample sex'] = df.ix[i, 'Biosample sex'] \
            if not pd.isnull(df.ix[i, 'Biosample sex']) else ""
        char['Biosample Age'] = df.ix[i, 'Biosample Age'] \
            if not pd.isnull(df.ix[i, 'Biosample Age']) else ""
        sample.characteristics = char
        human_encode_map[sample.id] = sample

    df = df[['File accession', 'Biosample term id', output_prefix.capitalize() + "_Description",
             'Controlled by', 'Input_Description',
             'Assay', 'Biosample term name', 'Biosample type', 'Biosample organism',
             'File download URL', 'Platform', 'Experiment target', 'Confidence']]

    df.columns = ['Data_ID', 'Study_ID', "Data_Description",
                  'Input', 'Input_Description',
                  'Sequencing_Protocol', 'Cell Line', 'Cell Type', 'Species',
                  'Raw Data', 'Instrument_Model', 'Experiment target/antibody', 'Confidence']
    df = df.set_index(['Data_ID'])

    return samples_df, df, human_encode_map

