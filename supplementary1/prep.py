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
    for cl in terms.columns:
        terms[cl] = terms[cl].str.lower()
    temp = {}
    for name, gr in terms.groupby('term1'):
        temp[name.capitalize()] = list(terms.loc[terms['term1'] == name, 'term2'].values)
    return temp


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

    except:
        print("Invalid URL")

