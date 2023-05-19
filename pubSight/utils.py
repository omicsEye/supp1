from Bio import Entrez
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import numpy as np
import time
import requests


def fetch(term, year, email, api_key, datetype='pdat'):
    """
    :param term: str, search query
    :param year: int, year of publication in YYYY format
    :param email: str, an email address
    :param datetype: Type of date used to limit a search. The allowed values vary between
    Entrez databases, but common values are 'mdat' (modification date), 'pdat' (publication date) and 'edat'
    (Entrez date). Generally an Entrez database will have only two allowed values for datetype.
    :return: count of pubication in pubmed
    """
    # Set the email address to identify yourself to the NCBI server
    if email:
        Entrez.email = email
    if api_key is not None:
        Entrez.api_key = api_key
        Entrez.sleep_between_tries = 5
    else:
        Entrez.sleep_between_tries = 15
    handle = Entrez.esearch(db="pubmed", term=term, retmax=0,
                            rettype='count', mindate=year,
                            maxdate=year, datetype=datetype)
    records = Entrez.read(handle)
    return int(records['Count'])


def read_terms(file_path, delimiter='\t'):
    df = pd.read_table(file_path, delimiter=delimiter, header=0, encoding='unicode_escape')
    print (df)
    return df


def get_year_data(search_query, year, email, api_key, datetype='pdat'):
    years = []
    count = []
    cn = 0
    try:
        while cn < 3:
            pub_count = fetch(term=search_query, year=year, email=email, api_key=api_key, datetype=datetype)
            years.append(year)
            count.append(int(pub_count))

            if int(pub_count) == 0:
                cn += 1
            else:
                cn = 0
            year -= 1
    except requests.exceptions.HTTPError as errh:
        print("Http Error:", errh)
    except requests.exceptions.ConnectionError as errc:
        print("Error Connecting:", errc)
    except requests.exceptions.Timeout as errt:
        print("Timeout Error:", errt)
    except requests.exceptions.RequestException as err:
        print("OOps: Something Else", err)
    df = pd.DataFrame(zip(years, count), columns=['year', 'count'])
    return df.loc[df.loc[:, 'count'] > 0, :]


def get_from_pd(data, year, email, api_key, datetype='pdat', write=False, report_dir='.'):
    for i in range(data.shape[0]):
        tmp = get_year_data(search_query=data.iloc[i, 2], year=year, email=email, api_key=api_key, datetype=datetype)
        tmp.loc[:, 'main_term'] = data.iloc[i, 0]
        tmp.loc[:, 'sub_term'] = data.iloc[i, 1]
        if i == 0:
            df_main = tmp.copy()
        else:
            df_main = pd.concat([df_main, tmp], ignore_index=True)
        print(data.iloc[i, 2], 'done')
        time.sleep(2)
    if write:
        df_main.to_csv(report_dir+'/pubmed_data.tsv', sep='\t', index=False)
    return df_main


def n_colors(n, colormap='cividis', custom_palette=None):
    """
    Utility for defining N evenly spaced colors across a color map or custom palette.
    :param n: number of colors to generate
    :param colormap: name of colormap to use (default: 'cividis')
    :param custom_palette: custom list of colors to use (default: None)
    :return: list of N colors
    """
    if custom_palette is not None:
        return custom_palette[:n]
    else:
        cmap = plt.get_cmap(colormap)
        cmap_max = cmap.N
        return [cmap(int(k * cmap_max / (n - 1))) for k in range(n)]


def pubmed_plot(data, colormap='cividis', custom_palette=None, group_legend=True, report_dir='.'):
    if group_legend:
        legend = False
    else:
        legend = True

    num_plots = len(set(data.loc[:, 'main_term']))
    num_cols = min(num_plots, 3)
    num_rows = (num_plots - 1) // num_cols + 1
    name_list = list(set(data.loc[:, 'main_term']))
    name_list.sort()

    min_year = data.loc[:, 'year'].min()
    max_year = data.loc[:, 'year'].max()
    temp_0 = pd.DataFrame(range(min_year, max_year + 1), columns=['year'])
    tech_list = list(set(data.loc[:, 'sub_term']))
    tech_list.sort()
    colors = n_colors(len(tech_list), colormap=colormap, custom_palette=custom_palette)
    color_pal = {}
    for n, tech in enumerate(tech_list):
        color_pal[tech] = colors[n]

    patch_list = []
    for key in tech_list:
        data_key = mpatches.Patch(color=color_pal[key], label=key)
        patch_list.append(data_key)

    fig = plt.figure(figsize=(7.2, (1.25 * num_rows / 2)))
    gs = GridSpec(num_rows, num_cols, wspace=0.0, hspace=0.0)

    cn = 0
    for i in range(num_plots):

        temp = data.loc[data.loc[:, 'main_term'] == name_list[cn], :]
        temp = pd.pivot_table(data=temp,
                              index=['year'],
                              columns=['sub_term'],
                              values='count').reset_index()
        temp = temp_0.merge(temp, how='left')
        temp.fillna(0, inplace=True)

        row_idx = i // num_cols
        col_idx = i % num_cols
        ax = fig.add_subplot(gs[row_idx, col_idx])
        ax.tick_params(axis='both', which='major', labelsize=6)

        try:
            temp.plot(x='year', kind='bar', stacked=True,
                      color=color_pal, ax=ax, legend=legend)
        except:
            pass

        plt.xlabel("")
        ymin, ymax = ax.get_ylim()
        ax.set_yticks(np.round(np.linspace(ymin, ymax, 5), 0))
        ax.xaxis.set_tick_params(labelsize=6)
        start, end = ax.get_xlim()
        #ax.xaxis.set_ticks(np.arange(start, end, .25))
        diff_years = int(end -start)
        ax.set_xticks(np.round(np.linspace(start, end, diff_years*10), 0))
        ax.spines['top'].set_linewidth(0.1)
        ax.spines['left'].set_linewidth(0.5)
        ax.spines['right'].set_linewidth(0.1)
        ax.spines['bottom'].set_linewidth(0.5)
        if legend:
            ax.legend(loc='lower left', fontsize=5, ncol=1)

        if i % 3 == 1:
            ax.tick_params(axis="y", direction="in", pad=-15)
        if i % 3 == 2:
            ax.yaxis.tick_right()
        if i < num_cols:
            ax.yaxis.get_major_ticks()[0].label1.set_visible(False)
            ax.set_xticklabels([])
        else:
            ax.yaxis.get_major_ticks()[0].label1.set_visible(False)
            ax.set_xticks(ax.get_xticks()[::5])
            ax.xaxis.label.set_visible(False)
            if i >= num_cols * (num_rows - 1):
                ax.xaxis.get_major_ticks()[0].label1.set_visible(False)
                ax.tick_params(axis="x", direction="out", pad=1)

        ax.xaxis.set_tick_params(labelsize=6)
        ax.yaxis.set_tick_params(labelsize=6)

        ax.text(.5, .85, name_list[cn], transform=ax.transAxes, ha="center", weight='bold', size=7)
        cn += 1

    fig.text(0.5, -0.05, 'Year', ha='center', fontsize=7, weight='bold')
    fig.text(-0.01, 0.6, 'Number of publications', va='center', rotation='vertical', fontsize=7, weight='bold')
    plt.tight_layout(pad=0.05)
    if group_legend:
        if num_plots % 3 == 0:
            plt.legend(handles=patch_list, bbox_to_anchor=(1, -0.05),
                       ncol=4, prop={'size': 6}, bbox_transform=fig.transFigure)
        else:
            plt.legend(handles=patch_list, bbox_to_anchor=(.7, .5),
                       ncol=2, prop={'size': 6}, bbox_transform=fig.transFigure)

    fig.savefig(report_dir + "/pubmed_fig.pdf", dpi=600, bbox_inches="tight")
    fig.savefig(report_dir + "/pubmed_fig.png", dpi=600, bbox_inches="tight")
    return fig



