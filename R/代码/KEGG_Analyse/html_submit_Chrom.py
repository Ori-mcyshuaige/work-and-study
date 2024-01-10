# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 14:05:46 2021
@author: fanjingjing
"""
# pip install selenium
from selenium import webdriver
import time
import os.path
from bs4 import BeautifulSoup
import requests
# import sys
import re
import pandas as pd

def kegg_online(species,keggHitMin,path_int,path_out):
    ## chrome浏览器
    chrome_options = webdriver.ChromeOptions()
    chrome_options.add_argument('--headless')
    chrome_options.add_argument('--disable-gpu')
    browser = webdriver.Chrome(chrome_options=chrome_options)
    
    browser.get('https://www.kegg.jp/kegg/mapper/color.html')
    browser.find_element_by_id("s_map").clear()
    browser.find_element_by_id("s_map").send_keys(species)###------species
    browser.find_element_by_name("color_list").send_keys(path_int)#file-----path_int
    browser.find_element_by_xpath('.//*[@value="Exec"]').submit()
    time.sleep(10)# 保证浏览器响应成功后再进行下一步操作

    html_file = browser.page_source
    browser.quit()
    
    soup = BeautifulSoup(html_file)
    # print(soup.prettify())
    
    #出表
    mystr=re.split('\n\n', soup.text)
    mystr2=[]
    for s in mystr:
        if s.startswith(species):
            s2=re.sub("\xa0", "", s)
            mystr2.append(s2)
            # print ("删除空值后的输出如下:\n",mystr2)
            
    mystr3=[]
    for s in mystr:
        if s.startswith("\xa0\xa0"):
            s2=re.sub("\n", "||", s)
            s3=re.sub("\xa0", "", s2)
            mystr3.append(s3)
            # print ("删除空值后的输出如下:\n",mystr2)
    
    if(len(mystr2)==len(mystr3)):#将两个list合并后输出
        mydata=[]
        for i in mystr2:#数据框第一列提取
            # print(i)
            num=re.findall(r"\d+",i)
            mydata.append(num[len(num)-1])
        z = list(zip(mystr2,mystr3,mydata))
        z2=pd.DataFrame(z,columns=['Pathway','Compound','num'])
        z2.num=pd.to_numeric(z2.num)
        out=z2[z2.num>=keggHitMin]
        del out['num']
        # z2.to_excel(os.path.join(path_out,"KEGG_Pathway.xlsx"), index=False)

        # Create a Pandas Excel writer using XlsxWriter as the engine.
        writer = pd.ExcelWriter(os.path.join(path_out,"KEGG_Pathway.xlsx"), engine='xlsxwriter')
        out.to_excel(writer, sheet_name='pathway', index=False)#追加sheet
        worksheet = writer.sheets['pathway']# Get the xlsxwriter worksheet object.
        worksheet.set_column("A:B", 80) # 在这里更改宽度值
        writer.save()# Close the Pandas Excel writer and output the Excel file.
    else:
        print("Please check the list!")
        return
        
    #出图
    args=[]
    for item in soup.find_all("a"):
        # print(item.get("href"))
        # print(item.text)
        ii = item.get("href")
        if "args" in ii:
            # print(ii)
            args.append(ii)
    k = [str(r"https://www.kegg.jp")+i for i in args]##列表生成式
    
    pic=[]###从上面出的表out,hsa来选图
    for p in list(out.Pathway):
        pic.append(re.split(" ", p)[0])
    
    k2=[]
    for i in pic:
        for j in k:
            if i in j:
                k2.append(j)
                
    driver = webdriver.Chrome(chrome_options=chrome_options)
    for i in k2:
        img_url = i
        driver.get(img_url)#driver.page_source
        # 网页中图片的定位(也可用其他方式定位)
        # img = driver.find_element_by_xpath("//img[@id='pathwayimage']")
        try:#报错执行,网页下包含非pathwayimage的图片
            img = driver.find_element_by_id("pathwayimage")
            img_url=img.get_attribute("src")
            file_name = img_url.split("/")[len(img_url.split("/"))-1]
            r = requests.get(img_url)
            with open(os.path.join(path_out,"pic",file_name), 'wb') as f:
                f.write(r.content)
        except Exception as e:
            print("----------------------------")
            print(e)
            print(i)
            print("This isn't a pathway image!!")

        

        
if __name__ == "__main__":    
    os.chdir(r"D:\3.好用的脚本收集\python\KEGG_Analyse")
    CUR_PATH = r'D:\3.好用的脚本收集\python\KEGG_Analyse'
    species="hsa"
    keggHitMin = int(1)
    for i in os.listdir(os.path.join(CUR_PATH,"input")):
        path_int=os.path.join(CUR_PATH,"input",i)
        path_out=os.path.join(CUR_PATH,"output",re.split("\\.",i)[0])
        if not os.path.exists(os.path.join(path_out,"pic")):
            os.makedirs(os.path.join(path_out,"pic"))
        kegg_online(species = species,keggHitMin = keggHitMin,
                    path_int=path_int,path_out=path_out)
        
