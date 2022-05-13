from prep import *
from viz import *

new_dir = make_dir()
dl_dir = os.getcwd() + '\\' + new_dir + '\\'
terms = read_terms2('terms.csv')

downloader(terms=terms, download_dir=dl_dir)

ls = sorted(list(set(terms['term1'])))

vizz(name_list=ls, report_dir=dl_dir)