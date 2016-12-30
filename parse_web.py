#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from multiprocessing import Pool
import requests
import bs4
import json
import io
import argparse
import time
import re

#root_url = "https://www.zhihu.com/"
root_url = "http://www.juzimi.com"

def get_q_url(single_id):    # this function generate one URL
    return root_url+'/ju/'+str(single_id)

def get_URLs(num):     # this function generate a series of URLs and store them in a list
    urls = map(get_q_url, range(20000, 20000 +num))
    return urls

def retrieve_data(url):
    dataList = {}
    user_agent = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/36.0.1985.143 Safari/537.36'
    header = {"User-Agent": user_agent}
    res = requests.get(url, headers=header)  # response
    if res.status_code != 200:
        return {'noValue': 'noContent'}

    soup = bs4.BeautifulSoup(res.text, 'html.parser')
    #dataList["index"] = soup.title.string[4:7]     # you could select other contents from the target html and store it in dataList

    for meta in soup.select('meta'):
        if meta.get('name') == 'description':
            dataList["content"] = meta.get("content")

    #dataList['imgUrl'] = soup.find_all('img')[1]['src']  # Of course you could store image information with corresponding URL
    return dataList

if __name__ == '__main__':
    pool = Pool(2)
    all_Urls = get_URLs(4)
    ResList = pool.map(retrieve_data, all_Urls)     # visit url in all_Urls and get all contents from them , then store them in ResList 
    jsonData = json.dumps({'data': ResList}, ensure_ascii=False, encoding='utf8')
    with io.open('data1.txt', 'w', encoding = 'utf-8') as Out:
        Out.write(jsonData)
        #json.dump(jsonData, Out) 

