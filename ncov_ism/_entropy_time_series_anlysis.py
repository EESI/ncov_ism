import logging
import pickle
import pandas as pd
import numpy as np
import datetime
from ._analyzeism import time_subset
from ._pickism import entropy_analysis

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

def entropy_time_series_analysis(INPUT_FOLDER, OUTPUT_FOLDER, REFERENCE_ID):
    '''
    perform entropy analysis on sequences collected in different time period to identify covarying positions.
    '''
    data_df_all = pd.read_pickle('{}/data_df_with_correction.pkl'.format(INPUT_FOLDER))
    data_df_all['date'] = pd.to_datetime(data_df_all['date'])

    REFERENCE = (REFERENCE_ID, data_df_all[data_df_all['gisaid_epi_isl'] == REFERENCE_ID]['sequence'].iloc[0])

    seq_index = []
    index = 0
    for base in REFERENCE[1]:
        if base == '-':
            seq_index.append(index)
        else:
            index += 1
            seq_index.append(index)

    reference_local_index_map = np.array(seq_index)
    MAX_DATE = data_df_all['date'].max().date()

    FLAG = True
    start = datetime.datetime.strptime('2019-12-24', "%Y-%m-%d")
    time_list = []
    while FLAG:
        end = start + datetime.timedelta(days=14)
        time_list.append((str(start.date()), str(end.date())))
        start = end
        if end.date() > MAX_DATE:
            FLAG = False

    for start, end in time_list:
        data_df = time_subset(data_df_all, '2019-12-24', end)
        logging.info('Entropy Time Series Analysis: processing data between {} and {} ({} sequences).'.format('2019-12-24', end, data_df.shape[0]))
        H_list, null_freq_list = entropy_analysis(data_df)
        H_list = np.array(H_list)
        null_freq_list = np.array(null_freq_list)
        pickle.dump([H_list, null_freq_list], open('{}/ENTS_{}_{}.pkl'.format(OUTPUT_FOLDER, '2019-12-24', end), 'wb'))

    logging.info('Entropy Time Series Analysis: DONE.')
