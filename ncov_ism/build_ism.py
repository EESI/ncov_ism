from ._loaddata import load_data
from ._pickism import entropy_analysis, pick_ISM_spots, annotate_ISM, ISM_disambiguation
from ._analyzeism import ISM_analysis

import os
import pandas as pd
import logging
logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

def build_from_existing(REFERENCE, H_list, OUTPUT_FOLDER):
    if not os.path.isfile('{}/ISM_annotation.txt'.format(OUTPUT_FOLDER)):
        return None
    
    annotation_df = pd.read_csv('{}/ISM_annotation.txt'.format(OUTPUT_FOLDER))
    
    ref_to_align = {}
    index = 0
    for idx, base in enumerate(REFERENCE[1]):
        if base != '-':
            ref_to_align[index] = idx
            index += 1
    
    min_en = 100
    for ref_idx in annotation_df['Ref position'].tolist():
        tmp_en = H_list[ref_to_align[ref_idx-1]]
        if tmp_en < min_en:
            min_en = tmp_en
            
    return min_en

def build_ISM(MSA_FILE_NAME, META_FILE_NAME, reference_genbank_name, OUTPUT_FOLDER, REFERENCE_ID, en_thres, null_thres):
    '''
    build ISMs from multiple sequence alignment and metadata file
    '''    
    logging.info('loading genomic sequences and metadata: ...')
    
    data_df = load_data(MSA_FILE_NAME, META_FILE_NAME)
    
    REFERENCE = (REFERENCE_ID, data_df[data_df['gisaid_epi_isl'] == REFERENCE_ID]['sequence'].iloc[0])
    REFERENCE_date = data_df[data_df['gisaid_epi_isl'] == REFERENCE_ID]['date'].min().date()

    H_list, null_freq_list = entropy_analysis(data_df)
    
    ## ===== choose entropy such that all old positions are preserved ===== ##
    en_min = build_from_existing(REFERENCE, H_list, OUTPUT_FOLDER)
    if en_min is not None and en_min < en_thres:
        en_thres = en_min - 0.0001
        logging.info('Informative Subtype Marker picking: using entropy threshold = {} to preserve all existing ISM positions'.format(en_thres))
    
    ## ===== choose entropy such that all old positions are preserved ===== ##
    
    position_list = pick_ISM_spots(H_list, null_freq_list, en_thres, null_thres)

    annotation_df = annotate_ISM(data_df, REFERENCE, position_list, reference_genbank_name)
    
    annotation_df.to_csv('{}/ISM_annotation.txt'.format(OUTPUT_FOLDER), index=False)
    
    data_df['ISM'] = data_df.apply(lambda x, position_list=position_list: ''.join([x['sequence'][position[0]] for position in position_list]), axis=1)
    ISM_df = data_df.drop(['sequence'], axis=1)

    data_df.to_pickle('{}/data_df_without_correction.pkl'.format(OUTPUT_FOLDER))
    ISM_df.to_csv('{}/ISM_df_without_correction.csv'.format(OUTPUT_FOLDER), index=False)

    logging.info('Informative Subtype Marker picking: ISM Table saved.')
    
    ISM_error_correction_partial, ISM_error_correction_full = ISM_disambiguation(ISM_df, THRESHOLD=0)

    ISM_df['ISM'] = ISM_df.apply(lambda x, 
             error_correction=ISM_error_correction_partial: error_correction[x['ISM']] if x['ISM'] in error_correction else x['ISM'],
             axis = 1)
    data_df['ISM'] = data_df.apply(lambda x, 
                 error_correction=ISM_error_correction_partial: error_correction[x['ISM']] if x['ISM'] in error_correction else x['ISM'],
                 axis = 1)

    data_df.to_pickle('{}/data_df_with_correction.pkl'.format(OUTPUT_FOLDER))
    ISM_df.to_csv('{}/ISM_df_with_correction.csv'.format(OUTPUT_FOLDER), index=False)
    acknowledgement_table = ISM_df[['gisaid_epi_isl', 'date', 'segment', 'originating_lab', 'submitting_lab', 'authors', 'url', 'date_submitted']]
    acknowledgement_table.to_csv('{}/acknowledgement_table.txt'.format(OUTPUT_FOLDER), index = False)
    return ISM_df
