######write to a txt
import os
from multiprocessing import Pool
import time
import pickle


def create_id2pw(onespecies):
    pw2id = {}
    print("species:", onespecies)
    parmeters_for_enrich = {}
    id2pw = open(r"/home/wang_xue/kegg_no_enrich/from_kegg_new/%s/elementid2pw.txt" % onespecies, "w", encoding = "utf-8")
    for onefile in os.listdir(r"/home/wang_xue/kegg_no_enrich/from_kegg_new/%s/maps/" % onespecies):
        if not onefile.endswith(".png"):
            with open(r"/home/wang_xue/kegg_no_enrich/from_kegg_new/%s/maps/%s"%(onespecies, onefile), encoding = "utf-8") as f:
                for oneline in f:
                    splitline = oneline.strip().split("\t")
                    try:
                        pw2id[onefile].add(splitline[0])
                    except:
                        pw2id[onefile] = set([splitline[0]])
                    if splitline[0].isdigit:
                        try:
                            parmeters_for_enrich[onefile].add(splitline[0])
                        except:
                            parmeters_for_enrich = set([splitline[0]])
    
    # write into elementid2pw.txt
    for onepw in pw2id:
        for oneid in pw2id[onepw]:
            id2pw.write(oneid + "\t" + onepw + "\n")
    id2pw.close()
    
    # for enrich analysis
    with open("/home/wang_xue/kegg_no_enrich/from_kegg_new/%s/parmeters_for_enrich.plk" % onespecies, "wb") as f:
        pickle.dump(parmeters_for_enrich, f)

    # for extract compound from elementid2pw.txt
    print("start", "   ", onespecies)
    f1 = open(r"/home/wang_xue/kegg_no_enrich/from_kegg_new/%s/compoundid2name.txt" % onespecies, "w", encoding = "utf-8")
    result = set()
    with open(r"/home/wang_xue/kegg_no_enrich/from_kegg_new/%s/elementid2pw.txt" % onespecies) as f:
        for oneline in f:
            splitline = oneline.split("\t")
            if oneline.startswith("C") and oneline[1:6].isdigit():
                compound = splitline[0]
                splitfirst = compound.split(" (")
                ID = splitfirst[0]
                NAME = splitfirst[1][:-1]
                result.add(ID + "\t" + NAME + "\n")
    for one_C in result:
        f1.write(one_C)
    f1.close()
    print("%s Done" % onespecies)

if __name__ == "__main__":
    s = time.time()
    all_species = os.listdir(r"/home/wang_xue/kegg_no_enrich/from_kegg_new/")
    for onespecies in all_species:
        create_id2pw(onespecies)
    print("done ", time.time() - s)
