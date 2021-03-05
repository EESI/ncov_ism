from ._loaddata import load_data
from ._pickism import entropy_analysis, pick_ISM_spots, annotate_ISM, ISM_disambiguation, ISM_disambiguate_fast
from ._analyzeism import ISM_analysis, time_subset

import os
import pickle
import pandas as pd
import numpy as np
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

def build_ISM_subset(MSA_FILE_NAME, META_FILE_NAME, start_date, end_date, reference_genbank_name, OUTPUT_FOLDER, REFERENCE_ID, en_thres, null_thres):
    '''
    build ISMs from multiple sequence alignment and metadata file
    '''
    logging.info('loading genomic sequences and metadata: ...')

    data_df = load_data(MSA_FILE_NAME, META_FILE_NAME)

    REFERENCE = (REFERENCE_ID, data_df[data_df['gisaid_epi_isl'] == REFERENCE_ID]['sequence'].iloc[0])
    REFERENCE_date = data_df[data_df['gisaid_epi_isl'] == REFERENCE_ID]['date'].min().date()

    data_df_subset = time_subset(data_df, start_date, end_date)
    logging.info('Number of sequences to be considered for new ISM hot spots: {}.'.format(data_df_subset.shape[0]))
    REF_to_POS = {}
    REF_IDX = 0
    for idx, base in enumerate(REFERENCE[1]):
        if base == '-':
            REF_to_POS[REF_IDX] = idx
        else:
            REF_IDX += 1
            REF_to_POS[REF_IDX] = idx

    H_list, null_freq_list = entropy_analysis(data_df_subset)
    new_position_list = pick_ISM_spots(H_list, null_freq_list, en_thres, null_thres)
    new_position_list = [pac[0] for pac in new_position_list]
    new_position_set = set(new_position_list)
    logging.info('Identified {} new ISM hot spots in the subset.'.format(len(new_position_set)))
    old_annotation = '{}/ISM_annotation.txt'.format(OUTPUT_FOLDER)
    try:
        annotation_df = pd.read_csv('{}/ISM_annotation.txt'.format(OUTPUT_FOLDER))
        old_ref_position_list = list(annotation_df['Ref position'])
        old_position_list = [REF_to_POS[REF_IDX] for REF_IDX in old_ref_position_list]
        new_position_set.update(old_position_list)
        logging.info('Found {} exsiting ISM hot spots in old annotation file.'.format(len(old_position_list)))
    except:
        logging.info('No existing ISM annotation files found for merging.')

    position_list = sorted(new_position_set)

    seq_index = []
    index = 0
    for base in REFERENCE[1]:
        if base == '-':
            seq_index.append(index)
        else:
            index += 1
            seq_index.append(index)
    reference_local_index_map = np.array(seq_index)
    ref_pos_list = {}
    for pos in position_list:
        if seq_index[pos] in ref_pos_list:
            continue
        ref_pos_list[seq_index[pos]] = pos
    position_list = sorted(ref_pos_list.values())

    H_list, null_freq_list = entropy_analysis(data_df)
    pickle.dump([H_list, null_freq_list], open('{}/overall_entropy.pkl'.format(OUTPUT_FOLDER), 'wb'))
    logging.info('Found {} overall ISM hot spots.'.format(len(position_list)))
    logging.info('Overall positions: {}.'.format(', '.join([str(idx) for idx in position_list])))
    position_list = [(base_idx, H_list[base_idx]) for base_idx in position_list]
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
    return ISM_df, [pac[0] for pac in position_list]

def build_ISM_subset_fast(MSA_FILE_NAME, META_FILE_NAME, start_date, end_date, reference_genbank_name, OUTPUT_FOLDER, REFERENCE_ID, en_thres, null_thres):
    '''
    build ISMs from multiple sequence alignment and metadata file
    '''
    logging.info('loading genomic sequences and metadata: ...')

    data_df = load_data(MSA_FILE_NAME, META_FILE_NAME)

    REFERENCE = (REFERENCE_ID, data_df[data_df['gisaid_epi_isl'] == REFERENCE_ID]['sequence'].iloc[0])
    REFERENCE_date = data_df[data_df['gisaid_epi_isl'] == REFERENCE_ID]['date'].min().date()

    data_df_subset = time_subset(data_df, start_date, end_date)
    logging.info('Number of sequences to be considered for new ISM hot spots: {}.'.format(data_df_subset.shape[0]))
    REF_to_POS = {}
    REF_IDX = 0
    for idx, base in enumerate(REFERENCE[1]):
        if base == '-':
            REF_to_POS[REF_IDX] = idx
        else:
            REF_IDX += 1
            REF_to_POS[REF_IDX] = idx

    H_list, null_freq_list = entropy_analysis(data_df_subset)
    new_position_list = pick_ISM_spots(H_list, null_freq_list, en_thres, null_thres)
    new_position_list = [pac[0] for pac in new_position_list]
    new_position_set = set(new_position_list)
    logging.info('Identified {} new ISM hot spots in the subset.'.format(len(new_position_set)))
    old_annotation = '{}/ISM_annotation.txt'.format(OUTPUT_FOLDER)
    try:
        annotation_df = pd.read_csv('{}/ISM_annotation.txt'.format(OUTPUT_FOLDER))
        old_ref_position_list = list(annotation_df['Ref position'])
        old_position_list = [REF_to_POS[REF_IDX] for REF_IDX in old_ref_position_list]
        new_position_set.update(old_position_list)
        logging.info('Found {} exsiting ISM hot spots in old annotation file.'.format(len(old_position_list)))
    except:
        logging.info('No existing ISM annotation files found for merging.')

    position_list = sorted(new_position_set)

    seq_index = []
    index = 0
    for base in REFERENCE[1]:
        if base == '-':
            seq_index.append(index)
        else:
            index += 1
            seq_index.append(index)
    reference_local_index_map = np.array(seq_index)
    ref_pos_list = {}
    for pos in position_list:
        if seq_index[pos] in ref_pos_list:
            continue
        ref_pos_list[seq_index[pos]] = pos
    position_list = sorted(ref_pos_list.values())

    H_list, null_freq_list = entropy_analysis(data_df)
    pickle.dump([H_list, null_freq_list], open('{}/overall_entropy.pkl'.format(OUTPUT_FOLDER), 'wb'))
    logging.info('Found {} overall ISM hot spots.'.format(len(position_list)))
    logging.info('Overall positions: {}.'.format(', '.join([str(idx) for idx in position_list])))
    position_list = [(base_idx, H_list[base_idx]) for base_idx in position_list]
    annotation_df = annotate_ISM(data_df, REFERENCE, position_list, reference_genbank_name)

    annotation_df.to_csv('{}/ISM_annotation.txt'.format(OUTPUT_FOLDER), index=False)

    data_df['ISM'] = data_df.apply(lambda x, position_list=position_list: ''.join([x['sequence'][position[0]] for position in position_list]), axis=1)
    ISM_df = data_df.drop(['sequence'], axis=1)

    data_df.to_pickle('{}/data_df_without_correction.pkl'.format(OUTPUT_FOLDER))
    ISM_df.to_csv('{}/ISM_df_without_correction.csv'.format(OUTPUT_FOLDER), index=False)

    logging.info('Informative Subtype Marker picking: ISM Table saved.')

    ISM_error_correction_partial = ISM_disambiguate_fast(set(ISM_df['ISM']))

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
    return ISM_df, [pac[0] for pac in position_list]
