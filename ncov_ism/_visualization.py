import logging
import matplotlib
import pickle
matplotlib.use('Agg')
import matplotlib.colors as mcolors
import numpy as np
import matplotlib.pyplot as plt
plt.ioff()
font = {# 'family' : 'serif', # Times (source: https://matplotlib.org/tutorials/introductory/customizing.html)
        'family': 'sans-serif', # Helvetica
        'size'   : 12}

matplotlib.rc('font', **font)
text = {'usetex': False}
matplotlib.rc('text', **text)
monospace_font = {'fontname':'monospace'}
CSS4_COLORS = mcolors.CSS4_COLORS

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

def ISM_filter(dict_freq, threshold):
    """
    collapse low frequency ISMs into "OTHER" per location
    Parameters
    ----------
    dict_freq: dictionary
        ISM frequency of a location of interest
    threshold: float
        ISMs lower than this threshold will be collapsed into "OTHER"
    Returns
    -------
    res_dict: dictionary
        filtered ISM frequency of a location of interest
    """
    res_dict = {'OTHER': [0, 0]}
    total = sum([int(dict_freq[ISM][1]) for ISM in dict_freq])
    for ISM in dict_freq:
        if int(dict_freq[ISM][1])/total < threshold:
            res_dict['OTHER'] = [0, res_dict['OTHER'][1] + int(dict_freq[ISM][1])]
        else:
            res_dict[ISM] = [dict_freq[ISM][0], int(dict_freq[ISM][1]) + res_dict.get(ISM, [0, 0])[1]]
    if res_dict['OTHER'][1] == 0:
        del res_dict['OTHER']
    return res_dict

def ISM_time_series_filter(dict_freq, threshold):
    """
    collapse low frequency ISMs into "OTHER" per location
    Parameters
    ----------
    dict_freq: dictionary
        ISM frequency of a location of interest
    threshold: float
        ISMs lower than this threshold will be collapsed into "OTHER"
    Returns
    -------
    res_dict: dictionary
        filtered ISM frequency of a location of interest
    """
    res_dict = {'OTHER': [0, 0]}
    total = sum([int(dict_freq[ISM]) for ISM in dict_freq])
    for ISM in dict_freq:
        if int(dict_freq[ISM])/total < threshold:
            res_dict['OTHER'] = [0, res_dict['OTHER'][1] + int(dict_freq[ISM])]
        else:
            res_dict[ISM] = [dict_freq[ISM], int(dict_freq[ISM]) + res_dict.get(ISM, [0, 0])[1]]
    if res_dict['OTHER'][1] == 0:
        del res_dict['OTHER']
    return res_dict

def ISM_visualization(region_raw_count, state_raw_count, count_dict, region_list, state_list, time_series_region_list,
                      output_folder, ISM_FILTER_THRESHOLD=0.05, ISM_TIME_SERIES_FILTER_THRESHOLD=0.025):
    '''
    Informative Subtype Marker analysis visualization
    Parameters
    ----------
    region_raw_count: dictionary
        ISM frequency per region
    state_raw_count: dictionary
        ISM frequency per state
    count_dict: dictionary
        ISM frequency time series per region
    region_list: list
        regions of interest
    state_list: list
        states of interest
    time_series_region_list: list
        regions of interest for time series analysis
    output_folder: str
        path to the output folder
    ISM_FILTER_THRESHOLD: float
        ISM filter threshold
    ISM_TIME_SERIES_FILTER_THRESHOLD: float
        ISM filter threshold for time series
    Returns
    -------
    Objects for downstream visualization
    '''
    ISM_set = set([])

    region_pie_chart = {}
    for idx, region in enumerate(region_list):
        dict_freq_filtered = ISM_filter(region_raw_count[region], ISM_FILTER_THRESHOLD)
        region_pie_chart[region] = dict_freq_filtered
        ISM_set.update(dict_freq_filtered.keys())

    state_pie_chart = {}
    for idx, state in enumerate(state_list):
        dict_freq_filtered = ISM_filter(state_raw_count[state], ISM_FILTER_THRESHOLD)
        state_pie_chart[state] = dict_freq_filtered
        ISM_set.update(dict_freq_filtered.keys())

    count_list = []
    date_list = []
    
    sorted_date = sorted(count_dict.keys())

    for date in sorted_date:
        dict_freq = {}
        for region in time_series_region_list:
            regional_dict_freq = count_dict[date][region]
            dict_freq_filtered = ISM_time_series_filter(regional_dict_freq, ISM_TIME_SERIES_FILTER_THRESHOLD )
            ISM_set.update(list(dict_freq_filtered.keys()))
            dict_freq[region] = dict_freq_filtered
        count_list.append(dict_freq)
        date_list.append(date)
    return ISM_set, region_pie_chart, state_pie_chart, count_list, date_list

def customized_ISM_visualization(region_raw_count, count_dict, region_list, output_folder, 
                                 ISM_FILTER_THRESHOLD=0.05, ISM_TIME_SERIES_FILTER_THRESHOLD=0.025):
    '''
    Informative Subtype Marker analysis visualization
    Parameters
    ----------
    region_raw_count: dictionary
        ISM frequency per region
    state_raw_count: dictionary
        ISM frequency per state
    count_dict: dictionary
        ISM frequency time series per region
    region_list: list
        regions of interest
    state_list: list
        states of interest
    time_series_region_list: list
        regions of interest for time series analysis
    output_folder: str
        path to the output folder
    ISM_FILTER_THRESHOLD: float
        ISM filter threshold
    ISM_TIME_SERIES_FILTER_THRESHOLD: float
        ISM filter threshold for time series
    Returns
    -------
    Objects for downstream visualization
    '''
    ISM_set = set([])

    region_pie_chart = {}
    for idx, region in enumerate(region_list):
        dict_freq_filtered = ISM_filter(region_raw_count[region], ISM_FILTER_THRESHOLD)
        region_pie_chart[region] = dict_freq_filtered
        ISM_set.update(dict_freq_filtered.keys())

    count_list = []
    date_list = []
    
    sorted_date = sorted(count_dict.keys())

    for date in sorted_date:
        dict_freq = {}
        for region in region_list:
            regional_dict_freq = count_dict[date][region]
            dict_freq_filtered = ISM_time_series_filter(regional_dict_freq, ISM_TIME_SERIES_FILTER_THRESHOLD )
            ISM_set.update(list(dict_freq_filtered.keys()))
            dict_freq[region] = dict_freq_filtered
        count_list.append(dict_freq)
        date_list.append(date)
    return ISM_set, region_pie_chart, count_list, date_list

def get_color_names(CSS4_COLORS, num_colors):
    '''
    Prepare colors for each ISM.
    '''
    bad_colors = set(['seashell', 'linen', 'ivory', 'oldlace','floralwhite', 
                      'lightyellow', 'lightgoldenrodyellow', 'honeydew', 
                      'mintcream', 'azure', 'lightcyan', 'aliceblue', 
                      'ghostwhite', 'lavenderblush'
                     ])

    by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(color))),
                             name)
                            for name, color in CSS4_COLORS.items())
    names = [name for hsv, name in by_hsv][14:]
    prime_names = ['red', 'orange', 'green', 'blue', 'gold', 
                 'lightskyblue', 'brown', 'black', 'pink',
                 'yellow']
    OTHER = 'gray'
    name_list = [name for name in names if name not in prime_names and name != OTHER and name not in bad_colors]   
    if num_colors > len(name_list) - 10:
        logging.info('NOTE: Repetitive colors for different ISMs (inadequate distinctive colors)')
        name_list = name_list + ceil(num_colors/len(name_list)) * name_list
    if num_colors > len(prime_names):
        ind_list = np.linspace(0, len(name_list), num_colors - 10, dtype = int, endpoint=False).tolist()
        color_names = prime_names + [name_list[ind] for ind in ind_list]
    else:
        color_names = prime_names[:num_colors]
    return color_names

def global_color_map(COLOR_DICT, ISM_list, out_dir): 
    '''
    Plot color-ISM map for reference.
    Adapted from https://matplotlib.org/3.1.0/gallery/color/named_colors.html
    '''
    ncols = 3
    n = len(COLOR_DICT)
    nrows = n // ncols + int(n % ncols > 0)

    cell_width = 1300
    cell_height = 100
    swatch_width = 180
    margin = 30
    topmargin = 40

    width = cell_width * 3 + 2 * margin
    height = cell_height * nrows + margin + topmargin
    dpi = 300

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height,
                        (width-margin)/width, (height-topmargin)/height)

    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()
    # ax.set_title(title, fontsize=24, loc="left", pad=10)
    ISM_list.append('OTHER')
    for i, name in enumerate(ISM_list):
        row = i % nrows
        col = i // nrows
        y = row * cell_height

        swatch_start_x = cell_width * col
        swatch_end_x = cell_width * col + swatch_width
        text_pos_x = cell_width * col + swatch_width + 50

        ax.text(text_pos_x, y, name, fontsize=14,
                fontname='monospace',
                horizontalalignment='left',
                verticalalignment='center')

        ax.hlines(y, swatch_start_x, swatch_end_x,
                  color=COLOR_DICT[name], linewidth=18)
    plt.savefig('{}/COLOR_MAP.png'.format(out_dir), bbox_inches='tight', dpi=dpi)
    plt.close(fig)
    
def func(pct, allvals):
    '''
    covert to absolute value for pie chart plot.
    '''
    absolute = int(round(pct/100.*np.sum(allvals)))
    return "{:d}".format(absolute)

def plot_pie_chart(sizes, labels, colors, ax):
    '''
    plot pie chart
    Adapted from https://matplotlib.org/3.1.1/gallery/pie_and_polar_charts/pie_and_donut_labels.html#sphx-glr-gallery-pie-and-polar-charts-pie-and-donut-labels-py
    '''
    wedges, texts, autotexts = ax.pie(sizes, autopct=lambda pct: func(pct, sizes), colors = colors, textprops=dict(color="w"))
    time_labels = ['-' if label == 'OTHER' else label.split(' ')[1] for label in labels]
    ax.legend(wedges, time_labels,
#           title="Oligotypes",
          loc="lower left",
          bbox_to_anchor=(0.8, 0, 0.5, 1))
    ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    return wedges, labels

def regional_growth_plot(region, ISM_df, REFERENCE_date, count_list, date_list, COLOR_DICT, OUTPUT_FOLDER):
    '''
    time series plot for a region of interest
    '''
    xlim_len = (ISM_df[ISM_df['country/region'] == region]['date'].max().date() - REFERENCE_date).days
    
    fig = plt.figure(figsize = (30, 15))  
    n = 4
    ax=plt.subplot(1, 1, 1)
    
    regional_total = []
    ISM_regional_set = set([])
    for i in range(len(count_list)):
        regional_dict_freq = count_list[i][region]
        regional_total.append(sum([regional_dict_freq[ISM][1] for ISM in regional_dict_freq]))
        ISM_regional_set.update(regional_dict_freq.keys())
    ISM_regional_list = []
    for ISM in ISM_regional_set:
        if ISM != 'OTHER':
            ISM_regional_list.append(ISM)
    NONOTHER = len(ISM_regional_list)
    if 'OTHER' in ISM_regional_set:
        ISM_regional_list.append('OTHER')
    
    to_plot = []
    ISM_to_plot_list = []
    for ISM in ISM_regional_list:
        ISM_regional_growth = []
            
        for i in range(len(count_list)):
            regional_dict_freq = count_list[i][region]
            if ISM in regional_dict_freq and regional_dict_freq[ISM][1]!= 0:
                ISM_regional_growth.append(regional_dict_freq[ISM][1]/regional_total[i])
            else:
                if ISM == 'OTHER':
                    other_count = sum([regional_dict_freq[ISM][1] for ISM in regional_dict_freq if ISM not in ISM_regional_set])
                    if regional_total[i] != 0:
                        ISM_regional_growth.append(other_count/regional_total[i])
                    else:
                        ISM_regional_growth.append(0)
                else:
                    ISM_regional_growth.append(0)
        to_plot.append(ISM_regional_growth)
        ISM_to_plot_list.append(ISM)
        # ax.plot(ISM_regional_growth, color = COLOR_DICT[ISM], label = ISM, linewidth = 4, marker = 'o', markersize = 4)
    
    to_plot = np.row_stack(to_plot)
    y_stack = np.cumsum(to_plot, axis=0)
    x = np.arange(y_stack.shape[1])
    for i in range(y_stack.shape[0]):
        ISM = ISM_to_plot_list[i]
        if i == 0:
            ax.fill_between(x, 0, y_stack[i,:], facecolor=COLOR_DICT[ISM], alpha=.5, label = ISM)
        else:
            ax.fill_between(x, y_stack[i-1], y_stack[i,:], facecolor=COLOR_DICT[ISM], alpha=.5, label = ISM)
    
    major_ticks = np.arange(0, len(date_list), 5)
    minor_ticks = np.arange(0, len(date_list))
    major_label = []
    
    for i in major_ticks.tolist():
        major_label.append(str(date_list[i]))

    ax.set_xticks(minor_ticks, minor=True)
    ax.set_xticks(major_ticks)
    ax.set_xticklabels(major_label)
    plt.setp(ax.get_xticklabels(), rotation=90)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.legend(
              loc="lower left",
              bbox_to_anchor=(1, 0, 0.5, 1),
              prop={'family': monospace_font['fontname']})
    plt.xlim([-1, xlim_len])
    plt.ylabel('Relative abundance')
    ax.grid(which='minor', alpha=0.3, linestyle='--')
    ax.grid(which='major', alpha=0.8)
    plt.savefig('{}/3_ISM_growth_{}.png'.format(OUTPUT_FOLDER, region), bbox_inches='tight')
    plt.close(fig)  
    
def ISM_plot(ISM_df, ISM_set, region_list, region_pie_chart, state_list, state_pie_chart, REFERENCE_date, time_series_region_list, count_list, date_list, OUTPUT_FOLDER):
    '''
    Generate figures for ISM analysis.
    '''
    ISM_index = {}
    idx = 0
    for ISM, counts in ISM_df['ISM'].value_counts().items():
        ISM_index[ISM] = idx
        idx += 1

    logging.info('{} ISMs will show up in the visualizations'.format(len(ISM_set)))
    ISM_list = []
    for ISM in ISM_set:
        if ISM == 'OTHER':
            continue
        ISM_list.append((ISM, ISM_index[ISM]))
    ISM_list = sorted(ISM_list, key = lambda x: x[1])
    ISM_list = [item[0] for item in ISM_list]

    color_map = get_color_names(CSS4_COLORS, len(ISM_list))
    COLOR_DICT = {}
    for idx, ISM in enumerate(ISM_list):
        COLOR_DICT[ISM] = color_map[idx]
    COLOR_DICT['OTHER'] = 'gray'
    pickle.dump(COLOR_DICT, open('COLOR_DICT.pkl', 'wb'))
    global_color_map(COLOR_DICT, ISM_list, OUTPUT_FOLDER)
    
    DPI = 100
    fig = plt.figure(figsize=(25, 15))   

    wedges_list = []
    for idx, region in enumerate(region_list):
        dict_freq = region_pie_chart[region]
        total = sum([dict_freq[ISM][1] for ISM in dict_freq])
        labels = []
        sizes = []
        colors = []

        for ISM in dict_freq:
            if ISM == 'OTHER':
                continue
            labels.append('{}: {}'.format(ISM, dict_freq[ISM][0]))
            colors.append(COLOR_DICT[ISM])
            sizes.append(dict_freq[ISM][1])
        if 'OTHER' in dict_freq:
            labels.append('OTHER')
            colors.append(COLOR_DICT['OTHER'])
            sizes.append(dict_freq['OTHER'][1])

        ax=plt.subplot(5, 5, idx+1)
        wedges, labels = plot_pie_chart(sizes, labels, colors, ax)
        ax.set_title(region)
        wedges_list.append((wedges, labels))

    labels_handles = {}
    handles_OTHER = None
    for wedges, labels in wedges_list:
        for idx, label in enumerate(labels):
            label = label.split(':')[0]
            if label == 'OTHER':
                handles_OTHER = [wedges[idx], label]
                continue
            if label not in labels_handles:
                labels_handles[label] = wedges[idx]
    if handles_OTHER:
        handles_list = list(labels_handles.values()) + [handles_OTHER[0]]
        labels_list = list(labels_handles.keys()) + [handles_OTHER[1]]
        fig.legend(
          handles_list,
          labels_list,
          bbox_to_anchor=(0.82, 0.25),
          bbox_transform=plt.gcf().transFigure,
          ncol=5,
          prop={'family': monospace_font['fontname']}
        )
    else:
        fig.legend(
          labels_handles.values(),
          labels_handles.keys(),
          bbox_to_anchor=(0.82, 0.25),
          bbox_transform=plt.gcf().transFigure,
          ncol=5,
          prop={'family': monospace_font['fontname']}
        )
    plt.savefig('{}/1_regional_ISM.png'.format(OUTPUT_FOLDER), bbox_inches='tight', dpi=DPI, transparent=True)
    plt.close(fig)
    
    fig = plt.figure(figsize=(25, 20))   

    subplot_y = int(np.sqrt(len(state_list)))
    subplot_x = int(np.sqrt(len(state_list))) + 1

    if subplot_x * subplot_y < len(state_list):
        subplot_y = subplot_x

    wedges_list = []
    for idx, state in enumerate(state_list):
        dict_freq = state_pie_chart[state]
        total = sum([dict_freq[ISM][1] for ISM in dict_freq])
        labels = []
        sizes = []
        colors = []

        for ISM in dict_freq:
            if ISM == 'OTHER':
                continue
            labels.append('{}: {}'.format(ISM, dict_freq[ISM][0]))
            colors.append(COLOR_DICT[ISM])
            sizes.append(dict_freq[ISM][1])
        if 'OTHER' in dict_freq:
            labels.append('OTHER')
            colors.append(COLOR_DICT['OTHER'])
            sizes.append(dict_freq['OTHER'][1])

        ax=plt.subplot(subplot_x, subplot_y, idx+1)
        wedges, labels = plot_pie_chart(sizes, labels, colors, ax)
        ax.set_title(state)
        wedges_list.append((wedges, labels))

    labels_handles = {}
    handles_OTHER = None
    for wedges, labels in wedges_list:
        for idx, label in enumerate(labels):
            label = label.split(':')[0]
            if label == 'OTHER':
                handles_OTHER = [wedges[idx], label]
                continue
            if label not in labels_handles:
                labels_handles[label] = wedges[idx]
    if handles_OTHER:
        handles_list = list(labels_handles.values()) + [handles_OTHER[0]]
        labels_list = list(labels_handles.keys()) + [handles_OTHER[1]]
        fig.legend(
          handles_list,
          labels_list,
          bbox_to_anchor=(0.82, 0.25),
          bbox_transform=plt.gcf().transFigure,
          ncol=5,
          prop={'family': monospace_font['fontname']}
        )
    else:
        fig.legend(
          labels_handles.values(),
          labels_handles.keys(),
          bbox_to_anchor=(0.82, 0.25),
          bbox_transform=plt.gcf().transFigure,
          ncol=5,
          prop={'family': monospace_font['fontname']}
        )
    plt.savefig('{}/2_intra-US_ISM.png'.format(OUTPUT_FOLDER), bbox_inches='tight', dpi=DPI, transparent=True)
    plt.close(fig)
    
    font = {'family': 'sans-serif', # Helvetica
            'size'   : 25}

    matplotlib.rc('font', **font) 

    for region in time_series_region_list:
        regional_growth_plot(region, ISM_df, REFERENCE_date, count_list, date_list, COLOR_DICT, OUTPUT_FOLDER)

def customized_ISM_plot(ISM_df, ISM_set, region_list, region_pie_chart, REFERENCE_date, count_list, date_list, OUTPUT_FOLDER):
    '''
    Generate figures for ISM analysis.
    '''
    ISM_index = {}
    idx = 0
    for ISM, counts in ISM_df['ISM'].value_counts().items():
        ISM_index[ISM] = idx
        idx += 1

    logging.info('{} ISMs will show up in the visualizations'.format(len(ISM_set)))
    ISM_list = []
    for ISM in ISM_set:
        if ISM == 'OTHER':
            continue
        ISM_list.append((ISM, ISM_index[ISM]))
    ISM_list = sorted(ISM_list, key = lambda x: x[1])
    ISM_list = [item[0] for item in ISM_list]

    color_map = get_color_names(CSS4_COLORS, len(ISM_list))
    COLOR_DICT = {}
    for idx, ISM in enumerate(ISM_list):
        COLOR_DICT[ISM] = color_map[idx]
    COLOR_DICT['OTHER'] = 'gray'
    pickle.dump(COLOR_DICT, open('COLOR_DICT.pkl', 'wb'))
    global_color_map(COLOR_DICT, ISM_list, OUTPUT_FOLDER)
    
    DPI = 100
    fig = plt.figure(figsize=(25, 15))   

    wedges_list = []
    for idx, region in enumerate(region_list):
        dict_freq = region_pie_chart[region]
        total = sum([dict_freq[ISM][1] for ISM in dict_freq])
        labels = []
        sizes = []
        colors = []

        for ISM in dict_freq:
            if ISM == 'OTHER':
                continue
            labels.append('{}: {}'.format(ISM, dict_freq[ISM][0]))
            colors.append(COLOR_DICT[ISM])
            sizes.append(dict_freq[ISM][1])
        if 'OTHER' in dict_freq:
            labels.append('OTHER')
            colors.append(COLOR_DICT['OTHER'])
            sizes.append(dict_freq['OTHER'][1])

        ax=plt.subplot(5, 5, idx+1)
        wedges, labels = plot_pie_chart(sizes, labels, colors, ax)
        ax.set_title(region)
        wedges_list.append((wedges, labels))

    labels_handles = {}
    handles_OTHER = None
    for wedges, labels in wedges_list:
        for idx, label in enumerate(labels):
            label = label.split(':')[0]
            if label == 'OTHER':
                handles_OTHER = [wedges[idx], label]
                continue
            if label not in labels_handles:
                labels_handles[label] = wedges[idx]
    if handles_OTHER:
        handles_list = list(labels_handles.values()) + [handles_OTHER[0]]
        labels_list = list(labels_handles.keys()) + [handles_OTHER[1]]
        fig.legend(
          handles_list,
          labels_list,
          bbox_to_anchor=(0.82, 0.25),
          bbox_transform=plt.gcf().transFigure,
          ncol=5,
          prop={'family': monospace_font['fontname']}
        )
    else:
        fig.legend(
          labels_handles.values(),
          labels_handles.keys(),
          bbox_to_anchor=(0.82, 0.25),
          bbox_transform=plt.gcf().transFigure,
          ncol=5,
          prop={'family': monospace_font['fontname']}
        )
    plt.savefig('{}/1_regional_ISM.png'.format(OUTPUT_FOLDER), bbox_inches='tight', dpi=DPI, transparent=True)
    plt.close(fig)
    
    font = {'family': 'sans-serif', # Helvetica
            'size'   : 25}

    matplotlib.rc('font', **font) 

    for region in region_list:
        regional_growth_plot(region, ISM_df, REFERENCE_date, count_list, date_list, COLOR_DICT, OUTPUT_FOLDER)
