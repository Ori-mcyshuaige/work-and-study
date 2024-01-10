#!/home/share/anaconda3/bin/python
# coding:utf-8

import pandas as pd
import xlsxwriter
from urllib.request import urlopen
from urllib.request import urlretrieve
import re
import time
import os
from random import randint
from PIL import Image
import shutil
import random
import sys 


class KEGG_Analysis:
    '''help kegg analysais by local database'''
    
    def __init__(self, input_path, output_dir, species, limit, local_kegg_path = "/home/share/kegg_no_enrich"):
        self.input_path = input_path
        self.output_path = output_dir
        self.species = species
        self.output_dir = output_dir
        self.limit = limit - 1
        self.local_kegg_path = local_kegg_path

    def create_dir(self):
        '''for other input files except those in report'''
        try:
            os.mkdir(self.output_path)
        except Exception as err:
            pass

    def map_pw2geneid_pw2gene(self):
        '''search our input gene one pathway picture and output two dictory pw2gene and pw2geneud'''
        gene2color = {}
        with open(self.input_path) as f:
            gene2color = {i.strip().split("\t")[0] : [i.strip().split("\t")[1]]  for i in f.readlines() if i.rstrip().split("\t")[0] != ""}
        pw2geneid = {} # for picture
        pw2gene = {}   # for excel output
        for one_gene in gene2color:
            with open(r"%s/from_kegg_new/%s/elementid2pw.txt" % (self.local_kegg_path, self.species)) as pws:
                for one_line in pws:
                    geneid_and_name = [] # for one line geneid (genename) 
                    split_line = one_line.strip().split("\t")
                    if " (" in split_line[0]:
                        geneid_and_name.append(split_line[0].split(" (")[0].strip())
                        geneid_and_name.append(split_line[0].split(" (")[1][:-1])
                    else:
                        geneid_and_name.append(split_line[0])
                        geneid_and_name.append(split_line[0])

                    if one_gene in geneid_and_name:
                        gene2color[one_gene].append(split_line[1])
                        # for download picture from kegg db
                        try:
                            pw2geneid[split_line[1]] += "/" + geneid_and_name[0] + "\t" + gene2color[one_gene][0] + ",black"
                        except:
                            pw2geneid[split_line[1]] = "/" + geneid_and_name[0] + "\t" + gene2color[one_gene][0] + ",black"
                        
                        # for create excel in report
                        try:
                            pw2gene[split_line[1]] += "||" + geneid_and_name[1]
                        except:
                            pw2gene[split_line[1]] = geneid_and_name[1]
        return(pw2gene, pw2geneid)
    
        
     
    def add_pw_description_and_gene_name(self, pw2gene):
        '''add pathway description and add color on significant gene or compound'''
        # create excel for report
        result_wb = xlsxwriter.Workbook(r"%s/KEGG Pathway.xlsx" % self.output_dir)
        sht = result_wb.add_worksheet("KEGG")
        sht.write(0, 0 , "KEGG pathway")
        sht.write(0, 1, "Genes")
        # for write into excel
        count = 2
        for one_pw in pw2gene:
            str1 = ""
            if pw2gene[one_pw].count("||") >= self.limit:
                genes = pw2gene[one_pw].split("||")
                count += 1
                # 1.for pathway description
                with open(r"%s/from_kegg_new/%s/pw2description.txt" % (self.local_kegg_path, self.species)) as pw2des:
                    for one_line in pw2des:
                        split_line = one_line.strip().split("\t")
                        if split_line[0][5:] == one_pw:
                            sht.write(count, 0, one_pw + ":" + split_line[1] + "(%s)" % (pw2gene[one_pw].count("||") + 1))
                            break
                # 2.for gene name and compound name
                n = 0 #for "||"
                for one_gene in genes:
                    # for gene name
                    with open(r"%s/from_kegg_new/%s/geneID2geneName.txt" % (self.local_kegg_path, self.species)) as geneID2geneName:
                        for one_line in geneID2geneName:
                            if ";" in one_line:
                                split_line = one_line.strip().split("\t")
                                symbol = split_line[1].split(";")[0].split(",")[0]
                                symbol_desc = split_line[1].split(";")[1]
                                gene_id = split_line[0]
                                if symbol == one_gene:
                                    n += 1
                                    if n >= 2:
                                        str1 += "||" + gene_id + " " + one_gene + ":" + symbol_desc
                                    else:
                                        str1 = gene_id + " " + one_gene + ":" + symbol + symbol_desc                           
                            # 更新的部分，有些输入的基因名为纯数字-------
                            else:
                                split_line = one_line.strip().split("\t")
                                symbolID = split_line[0].split(':')[1]
                                symbol = split_line[1]
                                symbol_desc = symbol
                                gene_id = split_line[0]
                                if symbolID == one_gene:
                                    n += 1
                                    if n >= 2:
                                        str1 += "||" + gene_id + " " + ":" + symbol_desc
                                    else:
                                        str1 = gene_id + " " + ":" + symbol + symbol_desc
                            
                    # for compound
                    with open(r"%s/from_kegg_new/%s/compoundid2name.txt" % (self.local_kegg_path, self.species)) as compound2Name:
                        for oneline in compound2Name:
                            splitline = oneline.strip().split("\t")
                            cpd_id = splitline[0]
                            cpd = splitline[1]
                            # id compound name = compound name
                            if cpd == one_gene:
                                n += 1
                                if n >= 2:
                                    str1 += "||" + "cpd:" + cpd_id + " " + cpd
                                else:
                                   str1 = "cpd: " + cpd_id + " " + cpd
                sht.write(count, 1, str1)
        result_wb.close()
        return(count)

        # download pw
    def copy_and_add_color(self, pw2geneid):
        for one_pw in pw2geneid:
            delete_logical = False
            if pw2geneid[one_pw].count("black") > self.limit:
                shutil.copy(r"%s/from_kegg_new/%s/maps/%s.png" % (self.local_kegg_path, self.species, one_pw), r"%s/%s.png" % (self.output_dir, one_pw))
                picture = Image.open(r"%s/%s.png" % (self.output_dir, one_pw))
                for one_element in pw2geneid[one_pw].split("/"):
                    if one_element != "":
                        with open(r"%s/from_kegg_new/%s/maps/%s"%(self.local_kegg_path, self.species, one_pw)) as coords:
                            for one_coord in coords:
                                geneid_in_coord = one_coord.strip().split("\t")[0].split(" ")[0].replace(":", "")
                                one_element_gene = one_element.split("\t")[0]
                                fourcoords = [int(x) for x in one_coord.strip().split("\t")[1].split(",")]
                                # color
                                showcolor = one_element.split("\t")[1].split(",")[0]
                                if one_element_gene == geneid_in_coord:
                                    if len(fourcoords) == 4: # for rect gene
                                        for x in range(fourcoords[0] + 1, fourcoords[2]):
                                            for y in range(fourcoords[1] + 1, fourcoords[3]):
                                                if picture.getpixel((x,y)) != (0,0,0,255):
                                                    if showcolor == "blue":
                                                        picture.putpixel((x, y), (0,0,255,255))
                                                    elif showcolor == "red":
                                                        picture.putpixel((x, y), (255, 0, 0, 255))
                                                    elif showcolor =="pink":
                                                        picture.putpixel((x, y), (255, 192, 203, 255))
                                    elif one_pw in [self.species + i for i in ["01100", "01110", "01120", "1130", "01200", "01210", "01212", "01230", "01220"]]: # for unfriendly pictures
                                        if os.path.isfile(r"%s/%s.png" % (self.output_dir, one_pw)):
                                            delete_logical = True
                                            break
                                    elif len(fourcoords) == 3:
                                        if showcolor == "blue":
                                            x1 = fourcoords[0]
                                            y1 = fourcoords[1]
                                            radius = fourcoords[2]
                                            for x in range(0, radius):
                                                for y in range(0, radius):
                                                    if picture.getpixel((x1 + x, y + y)) != (0, 0, 0, 255):
                                                        picture.putpixel((x1 + x, y1 + y), (0, 0, 255, 255)) 
                                            for x in range(0, radius):
                                                for y in range(0, radius):
                                                    if picture.getpixel((x1 - x, y1 - y)) != (0,0,0,255):
                                                        picture.putpixel((x1 - x, y1 - y), (0, 0, 255, 255))
                                            for x in range(0, radius):
                                                for y in range(0, radius):
                                                    if picture.getpixel((x1 + x, y1 - y)) != (0, 0, 0, 255):
                                                         picture.putpixel((x1 + x, y1 - y), (0,0,255,255))
                                            for x in range(0, radius):
                                                for y in range(0, radius):
                                                    if picture.getpixel((x1 - x, y1 + y)) != (0,0,0,255):
                                                        picture.putpixel((x1 - x, y1 + y), (0,0,255,255)) 
                                        elif showcolor == "red":
                                            x1 = fourcoords[0]
                                            y1 = fourcoords[1]
                                            radius = fourcoords[2]
                                            for x in range(0, radius):
                                                for y in range(0, radius):
                                                    if picture.getpixel((x1 + x, y + y)) != (0, 0, 0, 255):
                                                        picture.putpixel((x1 + x, y1 + y), (255, 0, 0, 255)) 
                                            for x in range(0, radius):
                                                for y in range(0, radius):
                                                    if picture.getpixel((x1 - x, y1 - y)) != (0,0,0,255):
                                                        picture.putpixel((x1 - x, y1 - y), (255, 0, 0, 255))
                                            # right and up
                                            for x in range(0, radius):
                                                for y in range(0, radius):
                                                     if picture.getpixel((x1 + x, y1 - y)) != (0, 0, 0, 255):
                                                         picture.putpixel((x1 + x, y1 - y), (255,0,0,255))
                                            # left and down
                                            for x in range(0, radius):
                                                for y in range(0, radius):
                                                    if picture.getpixel((x1 - x, y1 + y)) != (0,0,0,255):
                                                        picture.putpixel((x1 - x, y1 + y), (255,0,0,255)) 
                picture.save(r"%s/%s.png" % (self.output_dir, one_pw), "PNG")
                if delete_logical:
                    os.remove(r"%s/%s.png" % (self.output_dir, one_pw))
                
    def get_left_picture_from_kegg_new(self, count, pw2geneid):
        ''''download from kegg for some unfriendly picture'''
        ua_list = ["Mozilla/5.0 (Macintosh; Intel Mac OS X 10.6; rv2.0.1) Gecko/20100101 Firefox/4.0.1",
                "Mozilla/5.0 (Windows NT 6.1; rv2.0.1) Gecko/20100101 Firefox/4.0.1",
                "Opera/9.80 (Macintosh; Intel Mac OS X 10.6.8; U; en) Presto/2.8.131 Version/11.11",
                "Opera/9.80 (Windows NT 6.1; U; en) Presto/2.8.131 Version/11.11",
                "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_7_0) AppleWebKit/535.11 (KHTML, like Gecko) Chrome/17.0.963.56 Safari/535.11"]
        while len(os.listdir(self.output_dir)) < count + 1:
            for one_pw in pw2geneid:
                if (not os.path.isfile(r"%s/%s.png" % (self.output_dir, one_pw))) and (pw2geneid[one_pw].count("black") > self.limit):
                    regex = "/tmp/mark_pathway\d+?/%s\.png|/tmp/mark_pathway\d+?/%s_\d+?\.\d+?.png" % (one_pw, one_pw)
                    regexp = re.compile(regex)
                    my_headers = {'User-agent' : random.choice(ua_list)}
                    if not os.path.isfile(r"%s/%s.png" % (self.output_dir, one_pw)) and (pw2geneid[one_pw].count("black") > self.limit):
                        try:
                            #eg. https://www.kegg.jp/kegg-bin/show_pathway?map00400/1.14.16.1%09,blue/C00079%09,red/C00166/default%3dpink
                            img_url = r"https://www.kegg.jp/kegg-bin/show_pathway?%s%s" % (one_pw, pw2geneid[one_pw].replace("\t", "%09"))
                            html = urlopen(img_url)
                            html = str(html.read())
                            img_url = re.findall(regexp, html)[0]
                            urlretrieve("https://www.kegg.jp" + img_url, r"%s/%s.png" % (self.output_dir, one_pw))
                        except Exception as err:
                            print(err)
                            if os.path.isfile(r"%s/%s.png" % (self.output_dir, one_pw)):
                                os.remove(r"%s/%s.png" % (self.output_dir, one_pw))
                            continue
                    
if __name__ == "__main__":
    # for input
    local_path = os.path.abspath("~")
    input_path = "/home/share/kegg_no_enrich/input/"
    input_files = os.listdir(input_path)
    input_files_path = [input_path + i for i in input_files]
    # for output dirs
    output_dir = "/home/share/kegg_no_enrich/output/"
    output_dirs = [output_dir + i[:-4] for i in input_files]
    # run
    for i,j in enumerate(input_files_path):
        mykegg = KEGG_Analysis(input_path = j, output_dir = output_dirs[i], species = sys.argv[1], limit = 2)
        mykegg.create_dir()
        pw2gene,pw2geneid = mykegg.map_pw2geneid_pw2gene()
        count = mykegg.add_pw_description_and_gene_name(pw2gene)
        mykegg.copy_and_add_color(pw2geneid)
        mykegg.get_left_picture_from_kegg_new(count, pw2geneid)

