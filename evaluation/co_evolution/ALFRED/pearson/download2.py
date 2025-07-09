import requests
from bs4 import BeautifulSoup
import string
import os
import time

BASE_URL = "https://alfred.med.yale.edu/ALFRED/ALFREDServlet"
CHROMS = [str(i) for i in range(1,23)] + ["X", "PA"]
LETTERS = list(string.ascii_uppercase)
s = requests.Session()

def get_locus_uids(chrom, letter):
    s.post(BASE_URL, data={"alfred": "chromosome", "chrom": chrom})
    resp = s.post(BASE_URL, data={"alfred": "lociAlpha", "lociAlpha": letter})
    soup = BeautifulSoup(resp.text, "html.parser")
    table_div = soup.find("div", class_="CSSTableLong")
    if not table_div: return []
    table = table_div.find("table")
    uids = []
    for tr in table.find_all("tr")[1:]:
        tds = tr.find_all("td")
        if len(tds) < 1: continue
        form = tds[0].find("form")
        if form:
            locusUid = form.find("input", {"name": "locusUId"})["value"]
            uids.append(locusUid)
    return uids

def get_site_uids(locusUid, max_sites=6):
    resp = s.post(BASE_URL, data={"alfred": "locusInfo", "locusUId": locusUid})
    soup = BeautifulSoup(resp.text, "html.parser")
    site_table = soup.find("div", class_="CSSTableSite")
    if not site_table: return []
    table = site_table.find("table")
    siteUids = []
    for tr in table.find_all("tr")[1:]:
        tds = tr.find_all("td")
        form = tds[0].find("form")
        if form:
            siteUId = form.find("input", {"name": "siteUId"})["value"]
            siteUids.append(siteUId)
        if len(siteUids) >= max_sites:
            break
    return siteUids

def download_freq(siteUId, chrom):
    outdir = f"alfred_freqs2/{chrom}"
    os.makedirs(outdir, exist_ok=True)
    resp = s.post(BASE_URL, data={"alfred": "freqDownload", "siteUId": siteUId})
    path = f"{outdir}/{siteUId}_freq.txt"
    with open(path, "wb") as f:
        f.write(resp.content)
    print(f"Downloaded {path}")

for chrom in CHROMS:
    site_count = 0
    locus_checked = set()
    for letter in LETTERS:
        locus_uids = get_locus_uids(chrom, letter)
        for locusUid in locus_uids:
            if locusUid in locus_checked:
                continue
            locus_checked.add(locusUid)
            if site_count >= 100:
                break
            site_uids = get_site_uids(locusUid, max_sites=6)
            for siteUId in site_uids:
                if site_count >= 100:
                    break
                try:
                    download_freq(siteUId, chrom)
                    site_count += 1
                    time.sleep(0.2)
                except Exception as e:
                    print(f"Error {siteUId}: {e}")
        if site_count >= 100:
            break
    print(f"{chrom}: 共采集 {site_count} 个SNP位点频率文件。")