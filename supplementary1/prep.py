import os
import glob
import datetime
from shutil import move
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
import time
import pandas as pd
import urllib.parse as urlp


def make_dir():
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    new_dir = str('data_' + dt_label)
    os.makedirs(new_dir)
    return new_dir


def latest_file(dir_path, new_name):
    files = glob.glob(dir_path + r'/*')
    latest_file = max(files, key=os.path.getctime)
    new_name = dir_path+ '\\' + new_name + '.csv'
    move(latest_file, new_name)


def read_terms(file_name):
    terms = pd.read_csv(file_name)
    terms.fillna('', inplace=True)
    terms['term1'] = terms['term1'].str.capitalize()
    url_pre = 'https://pubmed.ncbi.nlm.nih.gov/?term='
    url_suf = '&sort=relevance'
    urls = []
    for i in range(len(terms)):
        tmp = terms.iloc[i, :].tolist()
        tmp = [x for x in tmp if len(x) > 1]
        tmp = " AND ".join(tmp)
        tmp = urlp.quote_plus(tmp)
        urls.append(url_pre + tmp + url_suf)
    terms['search_terms'] = urls
    return terms


def crawler(search_term, download_dir, sleep_time=1):
    # Headless option with window-size()
    chrome_options = webdriver.ChromeOptions()
    chrome_options.add_argument('--headless')
    chrome_options.add_argument('window-size=1920x1080')
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument('--disable-gpu')
    chrome_options.add_argument("--allow-insecure-localhost")
    chrome_options.add_experimental_option("prefs", {
        "download.default_directory": download_dir,
        "download.prompt_for_download": False
    })
    driver = webdriver.Chrome(chrome_options=chrome_options)
    try:
        driver.get(search_term)

        element = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.ID, "side-download-results-by-year-button"))
        )
        driver.execute_script("arguments[0].click();", element)
        # To verify that whether button is clicked or not
        driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")

        time.sleep(sleep_time)
        driver.close()

        check = True
    except:
        print("Invalid URL")
        check = False
    return check


def downloader(terms, download_dir):
    for i in range(len(terms)):
        search_term = terms.loc[i, 'search_terms']
        t1 = terms.loc[i, 'term1']
        t2 = terms.loc[i, 'term2']
        print('Downloading {} AN {}'.format(t1, t2))
        check = crawler(search_term=search_term, download_dir=download_dir)

        if check:
            latest_file(dir_path=download_dir, new_name=t1 + '_' + t2)
        else:
            print('There was no match for the terms')

    return print('Download Completed')
