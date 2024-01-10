
import sys
import pandas as pd
import xlsxwriter
from urllib.request import urlopen
from urllib.request import urlretrieve
from urllib.request import Request
import re
import time
import os
from random import randint
import random
from bs4 import BeautifulSoup
import pickle

with open("species.txt") as f:
    species_list = [i.strip() for i in f.readlines()]
# species_list = ["mmu"]
for one_species in species_list:
    print("one_species: ", one_species)
    location = r"/home/wang_xue/kegg_no_enrich"
    os.chdir(location)
    
    # download pathways from kegg
    if not os.path.exists(r"from_kegg_new/%s" % one_species):
        print("create %s dir" % one_species)
        os.chdir(r"%s/from_kegg_new" % location)
        os.mkdir(one_species)
        os.chdir(r"%s/from_kegg_new/%s" % (location, one_species))
    else:
        pass
    
    # download pathways list
    os.chdir(r"%s/from_kegg_new/%s" % (location, one_species))
    if not os.path.exists("maps"):
        os.mkdir("maps")
    if not os.path.isfile("pw2description.txt"):
        urlretrieve(r"http://rest.kegg.jp/list/pathway/%s" % one_species, "pw2description.txt")
    if not os.path.isfile("geneID2geneName.txt"):
        urlretrieve(r"http://rest.kegg.jp/list/%s" % one_species, "geneID2geneName.txt")
    
    
    os.chdir(location)
    
    
    
    # create all pathways list
    pathways = []
    with open("from_kegg_new/%s/pw2description.txt" % one_species) as f:
        for oneline in f:
            pathways.append(oneline.split("\t")[0][5:])
    
    
    ua_list = [
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.6; rv2.0.1) Gecko/20100101 Firefox/4.0.1",
            "Mozilla/5.0 (Windows NT 6.1; rv2.0.1) Gecko/20100101 Firefox/4.0.1",
            "Opera/9.80 (Macintosh; Intel Mac OS X 10.6.8; U; en) Presto/2.8.131 Version/11.11",
            "Opera/9.80 (Windows NT 6.1; U; en) Presto/2.8.131 Version/11.11",
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_7_0) AppleWebKit/535.11 (KHTML, like Gecko) Chrome/17.0.963.56 Safari/535.11"
    ]
    
    # delete empty file
    for oneEmptyFile in os.listdir(r"from_kegg_new/%s/maps" % one_species):
        print(oneEmptyFile)
        print(os.stat(r"from_kegg_new/%s/maps/%s" % (one_species, oneEmptyFile)).st_size)
        if os.stat(r"from_kegg_new/%s/maps/%s" % (one_species, oneEmptyFile)).st_size == 0:
            print(oneEmptyFile)
            os.remove(r"from_kegg_new/%s/maps/%s" % (one_species, oneEmptyFile))
    # download pictures
    count = len(pathways)
    
    while len(os.listdir(r"from_kegg_new/%s/maps" % one_species)) < 2 * count:
        for one_pw in pathways:
            pw2elementid = {} # for compound enrich
            print(one_pw)
            my_headers = {'User-agent' : random.choice(ua_list)}
            if not os.path.isfile(r"from_kegg_new\%s\maps\%s.png" % (one_species, one_pw)):
                print("Downloading ", one_pw)
                try:
                    print("picture %s starting!" % one_pw)
                    img_url = r"https://www.kegg.jp/kegg/pathway/%s/%s.png" % (one_species,one_pw)
                    #urlretrieve(img_url, "from_kegg_new/%s/maps/%s.png" % (one_species, one_pw))
                    img_request = Request(img_url, headers = my_headers)
                    img_response = urlopen(img_request, timeout = 60)
                    get_img = img_response.read()
                    with open("from_kegg_new/%s/maps/%s.png" % (one_species, one_pw), "wb") as fp:
                        fp.write(get_img)
                    print("picture %s done!" % one_pw)
                except Exception as err:
                        print(err)
                        if os.path.isfile("from_kegg_new/%s/maps/%s.png" % (one_species, one_pw)):
                            os.remove("from_kegg_new/%s/maps/%s.png" % (one_species, one_pw))
                        continue
            # loading geneName and compound name
            if not os.path.isfile(r"from_kegg_new/%s/maps/%s" % (one_species, one_pw)):
                print("file %s starting" % one_pw)
                try:
                    f = open(r"from_kegg_new/%s/maps/%s" % (one_species, one_pw), "w", encoding = "utf-8")
                    url1 = "https://www.kegg.jp/kegg-bin/show_pathway?%s" % one_pw
                    request = Request(url1, headers = my_headers)
                    response = urlopen(request, timeout = 60) # 超出时间主动放弃，进行下一个循环
                    bsObj = BeautifulSoup(response, "lxml")
                    for child in bsObj.findAll("area"):
                        for one_element in child.attrs["title"].split(", "):
                            #print(one_element)
                            f.write(one_element.strip() + "\t" + child.attrs["coords"].strip() + "\n")
                            try:
                                pw2elementid[one_pw].add(one_element.strip())
                            except:
                                pw2elementid[one_pw] = set(one_element.strip())
                    f.close()
                    print("file %s done" % one_pw)
                except Exception as err:
                    print(err)
                    f.close()
                    if os.path.isfile(r"from_kegg_new/%s/maps/%s" % (one_species, one_pw)):
                        os.remove(r"from_kegg_new/%s/maps/%s" % (one_species, one_pw))
                    continue
            for oneEmptyFile in os.listdir(r"from_kegg_new/%s/maps" % one_species):
                if os.stat(r"from_kegg_new/%s/maps/%s" % (one_species, oneEmptyFile)).st_size == 0:
                    os.remove(r"from_kegg_new/%s/maps/%s" % (one_species, oneEmptyFile))
            print("--------------------------------------------------------------------------")
