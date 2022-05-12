from prep import *
from viz import *

new_dir = make_dir()
dl_dir = os.getcwd() + '\\' + new_dir + '\\'
terms = read_terms('terms.csv')

for t1 in terms.keys():
    for t2 in terms[t1]:
        search_term = "https://pubmed.ncbi.nlm.nih.gov/?term=%28" + t1.replace(" ", "+") + \
              "%29+AND+%28%22" + t2.replace(" ", "+") + "%22%29&sort=relevance"
        print(search_term)

        crawler(search_term=search_term, download_dir=dl_dir)

        latest_file(dir_path=dl_dir, new_name=t1 + '_' + t2)

vizz(name_list=list(terms.keys()), report_dir=dl_dir)
