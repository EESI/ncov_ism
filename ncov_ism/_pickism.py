import logging
import json
import datetime
import pandas as pd
import numpy as np
import itertools
from collections import Counter, OrderedDict
from math import log2, ceil
from Bio import SeqIO

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

def base_entropy_masked(seq_list, base_set, base_idx):
    """
    Compute masked entropy for a position in the sequences
    Parameters
    ----------
    seq_list: list
        a list of aligned sequences
    base_set: set
        a set of unique characters
    base_idx: int
        a position in the sequences
    Returns
    -------
    H: float
        entropy
    masked_pct: float
        percentage of '-' and 'N'
    """
    # entropy analysis
    base_list = [seq[base_idx] for seq in seq_list]
    freq_dict = Counter(base_list)
    mask_list = ['-', 'N']
    n_seq = sum([freq_dict[base] for base in freq_dict if base not in mask_list])
    H = 0
    total_masked = 0
    for base in freq_dict:
        if base in mask_list:
            total_masked += freq_dict[base]
            continue
        P = freq_dict[base]/n_seq
        H -= log2(P) * P
    masked_pct = total_masked/len(base_list)
    return H, masked_pct

def entropy_analysis(data_df):
    """
    Masked Shannon entropy analysis for sequences
    Parameters
    ----------
    data_df: pandas.DataFrame
        merged Pandas dataframe
    Returns
    -------
    H_list: list
        entropy values for all positions
    null_freq_list: list
        masked percentage for all positions
    """
    seq_list = data_df['sequence'].values.tolist()
    base_set = set([])
    for seq in seq_list:
        base_set.update(set(seq))
    H_list = []
    null_freq_list = []
    STEP = ceil(len(seq_list[0]) / 10)
    for base_idx in range(len(seq_list[0])):
        if base_idx % STEP == 0:
            logging.info('Entropy analysis in progress: {}% completed.'.format(10 * base_idx // STEP))
        H, null_freq = base_entropy_masked(seq_list, base_set, base_idx)
        H_list.append(H,)
        null_freq_list.append(null_freq)
    logging.info('Entropy analysis in progress: DONE.')
    return H_list, null_freq_list

def entropy_analysis_customized(data_df, positions):
    """
    Masked Shannon entropy analysis for sequences
    Parameters
    ----------
    data_df: pandas.DataFrame
        merged Pandas dataframe
    Returns
    -------
    H_list: list
        entropy values for all positions
    null_freq_list: list
        masked percentage for all positions
    """
    seq_list = data_df['sequence'].values.tolist()
    base_set = set([])
    for seq in seq_list:
        base_set.update(set(seq))
    H_list = {}
    null_freq_list = {}
    STEP = ceil(len(positions) / 10)
    for base_idx in range(len(positions)):
        if base_idx % STEP == 0:
            logging.info('Entropy analysis in progress: {}% completed.'.format(10 * base_idx // STEP))
        H, null_freq = base_entropy_masked(seq_list, base_set, positions[base_idx])
        H_list[positions[base_idx]] = H
        null_freq_list[positions[base_idx]] = null_freq
    logging.info('Entropy analysis in progress: DONE.')
    return H_list, null_freq_list

def pick_ISM_spots(H_list, null_freq_list, en_thres=0.2, null_thres=0.25):
    """
    Pick Informative Subtype Markers based on masked Shannon entropy
    Parameters
    ----------
    H_list: list
        entropy values for all positions
    null_freq_list: list
        masked percentage for all positions
    en_thres: float
        Threshold for entropy values, entropy values greater than the threshold are valid
    null_thres: float
        Threshold for masked percentages, masked percentages lower than the threshold are valid
    Returns
    -------
    position_list: list
        selected position in a tuple: (postion, entropy values)
    """
    H_list, null_freq_list = np.array(H_list), np.array(null_freq_list)
    hot_spot_list = np.where((H_list > en_thres) & (null_freq_list < null_thres))[0]
    position_list = [(base_idx, H_list[base_idx]) for base_idx in hot_spot_list]
    logging.info('Pick Informative Subtype Markers: {} ISMs picked.'.format(len(position_list)))
    return position_list

def translate(seq): 
    '''
    Convert a given sequence of DNA into its Protein equivalent.
    Adapted from https://www.geeksforgeeks.org/dna-protein-python-3/
    Parameters
    ----------
    seq: str
        DNA sequence to translate
    Returns
    -------
    protein: str
        resultant prontein sequence 
    '''
    CODON_TABLE = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    
    protein = "" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= CODON_TABLE.get(codon, 'X')
    return protein 

def find_SNP(position, gene_dict, reference_raw):
    '''
    by Zhengqiao Zhao, v.0.2
    this function will take the position of a SNP as input,
    find the gene associated with this SNP if any
    and output the corresponding codon in reference. 
    Parameters
    ----------
    position: int
        1-indexed position of SNP
    gene_dict: dict
        a dictionary of gene sequences, key: (start, end), value: (name, amino acid sequence)
    reference_raw: str
        the nucleotide reference sequence
    Returns
    -------
    codon: str
        codon corresponding to the position
    codon_idx: int
        corresponding position in codon
    gene_name: str
        corresponding protein if valid
    '''
    for key in gene_dict:
    # interate over all genes and find the related gene
        if position > key[0] and position <= key[1] and (key[1] - key[0]) % 3 == 0:
            start = key[0]
            end = key[1]
            # extract the nucleotide gene sequence from the reference
            cDNA = reference_raw[start:end]
            # find the codon
            # python 0 indexed position
            delta = position -1 - key[0]
            condon_idx = delta % 3
            condon_pos = delta // 3 + 1
            full_codon = delta - condon_idx
            name, seq = gene_dict[(start, end)]
            
            if int(full_codon/3) >= len(seq):
                return None, None, name, None
            codon = cDNA[full_codon:full_codon+3]
            return codon, condon_idx, name, condon_pos
    return None, None, None, None

def load_gene_dict(reference_genbank_name="data/covid-19-genbank.gb"):
    """
    Load gene annotations from reference genbank file
    Parameters
    ----------
    reference_genbank_name: str
        path to the reference genbank file
    Returns
    -------
    gene_dict: dict
        dictionary containing gene annotations
    """
    recs = [rec for rec in SeqIO.parse(reference_genbank_name, "genbank")]
    gene_dict = {}
    for rec in recs:
        feats = [feat for feat in rec.features if feat.type == "CDS"]
        for feat in feats:
            content = '{}: {}'.format(feat.qualifiers['protein_id'][0], feat.qualifiers['product'][0])
            if feat.qualifiers['product'][0] == 'ORF1a polyprotein':
                continue
            if feat.location_operator == 'join':
                for item in feat.location.parts:
                    key = (item.start.position, item.end.position)
                    if 'translation' in feat.qualifiers:
                        seq = feat.qualifiers['translation']
                    if len(seq) == 1:
                        amino_acid_seq = seq[0]
                    gene_dict[key] = (content, amino_acid_seq)
            else:
                key = (feat.location.start.position, feat.location.end.position)
                if 'translation' in feat.qualifiers:
                    seq = feat.qualifiers['translation']
                if len(seq) == 1:
                    amino_acid_seq = seq[0]
                gene_dict[key] = (content, amino_acid_seq)
    return gene_dict

def annotate_ISM(data_df, REFERENCE, position_list, reference_genbank_name="data/covid-19-genbank.gb"):
    """
    Annotate Informative Subtype Markers
    Parameters
    ----------
    data_df: pandas.DataFrame
        merged Pandas dataframe
    REFERENCE: tuple
        tuple containing reference accession number and aligned reference genome
    position_list: list
        selected position in a tuple: (postion, entropy values)
    reference_genbank_name: str
        path to the reference genbank file
    Returns
    -------
    annotation_df: pandas.DataFrame
        Pandas dataframe containing annotaiton information
    """
    seq_list = data_df['sequence'].values.tolist()
    
    seq_index = []
    index = 0
    for base in REFERENCE[1]:
        if base == '-':
            seq_index.append(index)
        else:
            index += 1
            seq_index.append(index)
    reference_local_index_map = np.array(seq_index)
    mapped_reference_index = []
    for index, entropy in position_list:
        mapped_reference_index.append((index, reference_local_index_map[index], entropy))
    REFERENCE_ISM = ''.join([REFERENCE[1][item[0]] for item in position_list])
    logging.info('Reference ISM: {}.'.format(REFERENCE_ISM))
    
    gene_dict = load_gene_dict(reference_genbank_name)
    reference_raw = REFERENCE[1].replace('-', '')
    res = OrderedDict()
    res['Ref position'] = []
    res['Entropy'] = []
    res['Gene'] = []
    res['Is silent'] = []
    res['AA position'] = []
    for align_index, ref_index, entropy in mapped_reference_index:
        codon, codon_idx, name, codon_pos = find_SNP(ref_index, gene_dict, reference_raw)
        base_freq = Counter([item[align_index] for item in seq_list]).most_common()
        for alt_base, count in base_freq:
            if alt_base != reference_raw[ref_index-1]:
                break
        if codon is None:
            if_silence = True
        else:
            alt_codon = list(codon)
            alt_codon[codon_idx] = alt_base
            alt_codon = ''.join(alt_codon)
            ref_aa = translate(codon)
            ism_aa = translate(alt_codon)
            if ref_aa == ism_aa:
                if_silence = True
            else:
                if_silence = False
        res['Ref position'].append(ref_index)
        res['Entropy'].append(entropy)
        if name is None:
            name = 'Non-coding'
        res['Gene'].append(name)
        res['Is silent'].append(if_silence)
        if codon_pos is None:
            res['AA position'].append('NaN')
        else:
            res['AA position'].append('{}{}{}'.format(ref_aa, codon_pos, ism_aa))
    annotation_df = pd.DataFrame.from_dict(res)
    return annotation_df

# ambiguous bases correction
def is_same(error, target, mask, ISM_LEN):
    """
    Check if a masked ambiguous ISM is the same as a non-ambiguous ISM
    Parameters
    ----------
    error: str
        an ambiguous ISM
    target: str
        a non-ambiguous ISM
    mask: list
        positions masked
    ISM_LEN: int
        length of ISM
    Returns
    -------
    res: boolean
        if two ISMs are the same or not
    """
    match = np.array(list(target)) == np.array(list(error))
    res = np.logical_or(mask, match).sum() == ISM_LEN
    return res

def error_correction(error, ambiguous_base, base_to_ambiguous, ISM_list, ISM_LEN, THRESHOLD = 0.9):
    """
    Correct ISM by replacing ambiguous bases in an ISM by the similar non-ambiguous ISMs
    Parameters
    ----------
    error: str
        an ambiguous ISM
    ambiguous_base: set
        set containing  ambiguous bases
    base_to_ambiguous: dictionary
        map an ambiguous base to all possible bases
    ISM_list: list
        list containing all ISMs
    ISM_LEN: int
        length of ISM
    THRESHOLD: float
        percentage of non-ambiguous supporting instances
    Returns
    -------
    FLAG: boolean
        if fully corrected
    corrected: str
        corrected ISM
    """
    mask = [True if base in ambiguous_base else False for base in error]
    support_ISM = []
    for target_ISM in ISM_list:
        if is_same(error, target_ISM, mask, ISM_LEN):
            support_ISM.append(target_ISM)
    partial_correction = list(error)
    FLAG = True
    for position_idx in list(np.where(mask)[0]):
        possible_bases = set([candid_ISM[position_idx] for candid_ISM in support_ISM])
        possible_bases.discard('N')
        possible_bases.discard(error[position_idx])
        possible_bases.discard('-')
        non_ambiguous_set = set([])
        ambiguous_set = set([])
        for base in possible_bases:
            if base not in ambiguous_base:
                non_ambiguous_set.add(base)
            else:
                ambiguous_set.add(base)
        if len(ambiguous_set) == 0:
            if len(non_ambiguous_set) == 0:
                continue
            bases = ''.join(sorted(non_ambiguous_set))
            if len(bases) == 1:
                num_support = len([candid_ISM[position_idx] for candid_ISM in support_ISM if candid_ISM[position_idx] == bases])
                non_support = set([candid_ISM[position_idx] for candid_ISM in support_ISM if candid_ISM[position_idx] != bases])
                if num_support/len(support_ISM) > THRESHOLD and bases in ambiguous_base[error[position_idx]]:
                    partial_correction[position_idx] = bases
                else:
                    FLAG = False
                    logging.debug('Error Correction DEBUG: one-base-correction failed because no enough support: {}/{}: {}->{}'.format(num_support, len(support_ISM), non_support, bases))
            elif bases in base_to_ambiguous:
                FLAG = False
                partial_correction[position_idx] = base_to_ambiguous[bases]
            else:
                FLAG = False
                logging.debug("Error Correction DEBUG: can't find: {}".format(bases))
        else:
            bases_from_ambiguous_set = set([])
            ambiguous_bases_intersection = ambiguous_base[error[position_idx]].copy()  
            for base in ambiguous_set:
                bases_from_ambiguous_set = bases_from_ambiguous_set.union(ambiguous_base[base])
                ambiguous_bases_intersection = ambiguous_bases_intersection.intersection(ambiguous_base[base])
            
            if bases_from_ambiguous_set.issubset(ambiguous_base[error[position_idx]]) is False:
                logging.debug('Error Correction DEBUG: new bases {} conflict with or are not as good as original bases {}'.format(bases_from_ambiguous_set, ambiguous_base[error[position_idx]]))
                bases_from_ambiguous_set = ambiguous_base[error[position_idx]]
                
            bases_from_ambiguous_set = ''.join(sorted(bases_from_ambiguous_set))
            
            bases = ''.join(sorted(non_ambiguous_set))
            
            if len(bases) == 0:
                bases = bases_from_ambiguous_set
            if len(bases) == 1 and bases in bases_from_ambiguous_set:
                num_support = len([candid_ISM[position_idx] for candid_ISM in support_ISM if candid_ISM[position_idx] == bases])
                non_support = set([candid_ISM[position_idx] for candid_ISM in support_ISM if candid_ISM[position_idx] != bases])
                if num_support/len(support_ISM) > THRESHOLD and bases in ambiguous_bases_intersection:
                    partial_correction[position_idx] = bases
                else:
                    if bases not in ambiguous_bases_intersection:
                        logging.debug('Error Correction DEBUG: conflicts dected between proposed correct and all supporting ISMs')
                        bases = ''.join(ambiguous_bases_intersection.add(base))
                        if bases in base_to_ambiguous and set(bases).issubset(ambiguous_base[error[position_idx]]):
                            FLAG = False
                            partial_correction[position_idx] = base_to_ambiguous[bases]
                    else:
                        FLAG = False
                        logging.debug('Error Correction DEBUG: one-base-correction failed because no enough support: {}/{}: {}->{}'.format(num_support, len(support_ISM), non_support, bases))
            else:
                bases = ''.join(sorted(set(bases_from_ambiguous_set + bases)))
                if bases in base_to_ambiguous and set(bases).issubset(ambiguous_base[error[position_idx]]):
                    FLAG = False
                    partial_correction[position_idx] = base_to_ambiguous[bases]
                else:
                    FLAG = False
                    logging.debug('Error Correction DEBUG: new bases {} conflict with or are not as good as original bases {}'.format(bases, ambiguous_base[error[position_idx]]))
            
    return FLAG, ''.join(partial_correction)

def check_completeness(ISM):
    """
    Check if an ISM is fully corrected (no ambiguous bases)
    Parameters
    ----------
    ISM: str
        an ISM of interest
    Returns
    -------
    FLAG: boolean
        if fully corrected
    """
    for item in ISM:
        if item not in ['A', 'T', 'C', 'G', '-']:
            return False
    return True
    
def ISM_disambiguation(ISM_df, THRESHOLD=0):
    """
    Correct all ISMs by replacing ambiguous bases in an ISM by the similar non-ambiguous ISMs
    Ambiguous bases can be found in https://www.bioinformatics.nl/molbi/SCLResources/sequence_notation.htm
    Parameters
    ----------
    ISM_df: pandas.DataFrame
        Pandas dataframe containing an ISM column to be cleaned
    ambiguous_base: set
        set containing  ambiguous bases
    base_to_ambiguous: dictionary
        map an ambiguous base to all possible bases
    ISM_list: list
        list containing all ISMs
    ISM_LEN: int
        length of ISM
    THRESHOLD: float
        percentage of non-ambiguous supporting instances
    Returns
    -------
    FLAG: boolean
        if fully corrected
    corrected: str
        corrected ISM
    """
    
    ambiguous_base = {'B': set(['C', 'G', 'T']), 
                      'D': set(['A', 'G', 'T']),
                      'H': set(['A', 'C', 'T']),
                      'K': set(['G', 'T']),
                      'M': set(['A', 'C']),
                      'N': set(['A', 'C', 'G', 'T']),
                      'R': set(['A', 'G']),
                      'S': set(['C', 'G']),
                      'V': set(['A', 'C', 'G']),
                      'W': set(['A', 'T']),
                      'Y': set(['C', 'T'])}
    base_to_ambiguous = {}
    for base in ambiguous_base:
        bases = ''.join(sorted(ambiguous_base[base]))
        base_to_ambiguous[bases] = base

    ISM_list = list(ISM_df['ISM'].values)
    error_ISM_list = list(ISM_df[ISM_df.apply(lambda x, 
                              ambiguous_base=ambiguous_base: True if len(set(x['ISM']).intersection(ambiguous_base)) > 0 else False, 
                              axis = 1)]['ISM'].unique())
    ERR_DICT = {}
    for ISM in error_ISM_list:
        ERR_DICT[ISM] = ISM_df[ISM_df['ISM'] == ISM].shape[0]
    ISM_LEN = len(ISM_list[0])
    partial_ISM = 0
    partial_subj = 0
    full_ISM = 0
    full_subj = 0
    total_ISM = len(error_ISM_list)
    total_subj = sum([ERR_DICT[item] for item in ERR_DICT])
    ISM_error_correction_partial = {}
    ISM_error_correction_full = {}

    STEP = ceil(len(error_ISM_list) / 10)
    for ISM_idx, error in enumerate(error_ISM_list):
        if ISM_idx % STEP == 0:
            logging.info('ISM Disambiguation in progress: {}% completed.'.format(10 * ISM_idx // STEP))
        FLAG, correction = error_correction(error, ambiguous_base, base_to_ambiguous, ISM_list, ISM_LEN, THRESHOLD)
        FLAG = check_completeness(correction)
        if error != correction:
            ISM_error_correction_partial[error] = correction
            partial_ISM += 1
            partial_subj += ERR_DICT[error]
        if FLAG and error != correction:
            ISM_error_correction_full[error] = correction
            full_ISM += 1
            full_subj += ERR_DICT[error]
    logging.info('ISM Disambiguation in progress: DONE.')
    logging.info('ISM Disambiguation: percentage of unique ISMs partially corrected: {}'.format(partial_ISM/total_ISM))
    logging.info('ISM Disambiguation: percentage of unique ISMs completely corrected: {}'.format(full_ISM/total_ISM))
    logging.info('ISM Disambiguation: percentage of records (submissions) partially corrected: {}'.format(partial_subj/total_subj))
    logging.info('ISM Disambiguation: percentage of records (submissions) completely corrected: {}'.format(full_subj/total_subj))
    return ISM_error_correction_partial, ISM_error_correction_full


def ISM_disambiguate_fast(ism_list):
    '''Returns a dictionary mapping each distinct ISM to the corresponding disambiguated ISM, or if
    no disambiguation could be identified or it was already unambiguous, itself.
    If ismlist is a list nd not a set, then it instead returns a list of corrected ISMs,
    corresponding in order to original list of ISMs.
    '''
    if type(ism_list) is list:
        return_list = True
    else:
        return_list = False

    ism_set = set(ism_list)
    ismlen = len(next(iter(ism_list)))  # assume all ISMs are the same length


    ism_clean = set(); ism_toclean = set()
    # Create mask vectors for each ISM that is set to True if ACTG or gap vs. False if an ambiguous base,
    # and separate ISMs in the ones that ar e clean and the ones that are ambiguous.
    ism_mask_map = {}
    for ism in ism_set:
        ism_mask = [False]*ismlen
        for m in range(ismlen):
            if ism[m] in 'ACTG-':
                ism_mask[m] = True
        if any(ism_mask):
            ism_toclean.add(ism)
        else:
            ism_clean.add(ism)
        ism_mask_map[ism] = ism_mask

    ism_disamb = {}
    # assign all of the clean ISMs to themselves (i.e. they translate to themselves)
    # it is a one-element array since other ISMs may translate to multiple clean ISMs
    for ism in ism_clean:
        ism_disamb[ism] = ism


    def hamming_distance(string1, string2):
        # Return the Hamming distance between string1 and string2 of same length.
        distance = 0
        for i in range(len(string1)):
            if string1[i] != string2[i]:
                distance += 1
        # Return the final count of differences
        return distance


    # Find the Hamming distance of the ISM to each of the other ISMs
    # This means that *any* ISM, including one with ambiguous bases itself, may be a support
    HDist_Map = {}           # Dictionary where HDist_Map[ISM1][ISM2] = distance(ISM1,ISM2)
    for ism in ism_toclean:
        HDist_Map[ism] = {}
    # First, find all the binary distances between the ISMs that are to be cleaned
    for ismpair in itertools.combinations(ism_toclean, 2):
        HDist_Map[ismpair[0]][ismpair[1]] = hamming_distance(ismpair[0], ismpair[1])
        HDist_Map[ismpair[1]][ismpair[0]] = HDist_Map[ismpair[0]][ismpair[1]]
    # Second, find all the binary distances between each ISM to be cleaned with each clean ISM
    for ism in ism_toclean:
        for ismref in ism_clean:
            HDist_Map[ism][ismref] = hamming_distance(ism, ismref)

    for ism in ism_toclean:
        min_h_dist = min(HDist_Map[ism].values())
        min_hd = [key for key in HDist_Map[ism] if HDist_Map[ism][key]==min_h_dist]
        support = set(min_hd)
        # For each ambiguous base position in the ISM, see if there is a consensus among support set members
        # who have an ambiguity at that location.
        # Select the ambiguous base positions:
        amb_ntpos = set(itertools.compress(range(len(ism_mask_map[ism])), [not(m) for m in ism_mask_map[ism]]))
        # Identify the non-ambiguous bases in support ISMs at each of these positions. Once identified, overwrite
        # the position in a new disambiguated ISM:
        ism_new = ism
        for ntpos in amb_ntpos:
            supportbases = [supportism[ntpos] for supportism in support if ism_mask_map[supportism][ntpos]]
            if len(supportbases) == 0:
                # no valid bases at that position in any support ISMs
                pass
            elif supportbases.count(supportbases[0]) == len(supportbases):
                # then all of the bases in support ISMs are identical; replace the ambiguous base
                ism_new = ism_new[:ntpos] + supportbases[0] + ism_new[ntpos+1:]
        ism_disamb[ism] = ism_new

    if return_list:
        return [ism_disamb[ism] for ism in ism_list]
    else:
        return ism_disamb

