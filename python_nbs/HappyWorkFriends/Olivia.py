# -*- coding: utf-8 -*-
from .InviteFriends import *

def sigmoid(X):
    return( 1.0 / (1 + np.exp(-X)) )

def fpkm2tpm(fpkms):
    return(exp(log(fpkms) - log(sum(fpkms)) + log(1e+6)))

def FPKM2TPM(FPKMs):
    if not isinstance(FPKMs, (np.ndarray, np.generic)):
        FPKMs = np.array(FPKMs)
    return( FPKMs/sum(FPKMs) * 1e+6 )

def FPKM2TPM_log2unit(FPKMs):
    if not isinstance(FPKMs, (np.ndarray, np.generic)):
        FPKMs = np.array(FPKMs)
    FPKMs = np.exp2(FPKMs) - 1
    return( FPKMs/sum(FPKMs) * 1e+6 )

def count2tpm(counts, effLens):
    rate = log(counts) - log(effLens)
    denom = log(sum(exp(rate)))
    return(exp(rate - denom + log(1e6)))

def read_gene_effLen_and_hgncid(file='/mnt/afs/neo21c/neo/danxu/reference/GRCh37.p13/gencode.v22.annotation.gene.probeMap'):
    df = pd.read_table(file)
    df['effLen'] = df['chromEnd'] - df['chromStart']
    return( dict(zip(df['gene'], df['effLen'])), dict(zip(df['id'], df['gene'])) )

def write_tcga_gene_rnaexp(rnafile='/mnt/afs/neo21c/neo/danxu/reference/GRCh37.p13/TCGA-LUAD.htseq_fpkm.tsv', ofile='/home/danxu/2.files/database/tcga_gene_tpm.json'):
    df_rna = pd.read_table(rnafile)
    sample_cols = df_rna.columns[1:]
    dic_gene_effLen, dic_ensembl2hgnc = read_gene_effLen_and_hgncid()
    df_rna['Hugo_Symbol'] = [dic_ensembl2hgnc[ID] for ID in df_rna['Ensembl_ID']]
    #effLens = [dic_gene_effLen[gene] for gene in df_rna['Ensembl_ID']]
    df_rna[sample_cols] = df_rna[sample_cols].apply(FPKM2TPM_log2unit)
    dic_rnaexp_median = {}
    for n,i in df_rna.iterrows():
        dic_rnaexp_median[i['Hugo_Symbol']] = i[sample_cols].median().round(2)
    with open(ofile, 'w') as fo:
        json.dump(dic_rnaexp_median, fo, ensure_ascii=False)
        
def load_tcga_gene_rnaexp_dic(input_json_file='/home/danxu/2.files/database/tcga_gene_tpm.json'):
    with open(input_json_file, 'r') as f:
        return(json.load(f))

def read_and_reset_mhcpan_result_DFcolumns_base(file, start_cols=['Pos','Peptide','ID'], rep_cols=['core','icore','EL-score','EL_Rank','BA-score','BA_Rank']):
    df = pd.read_table(file, header=[0,1]).reset_index(drop=True)
    HLAs = open(file, 'r').readline().strip()
    ss = '\t' * len(rep_cols)
    HLAs = [i.replace('\t','') for i in HLAs.split(ss)]
    cols = [(i, "") for i in start_cols]
    for HLA in HLAs:
        for rep_col in rep_cols:
            cols.append((HLA, rep_col))
    cols.append(('Ave', '')); cols.append(('NB', ''))
    df.columns = pd.MultiIndex.from_tuples(cols)
    return(df)

def read_and_reset_mhcpan2_result_DFcolumns_base(file, start_cols=['Pos','Peptide','ID','Target'], rep_cols=['Score','Rank','Score_BA','nM','Rank_BA']):
    df = pd.read_table(file, header=[0,1]).reset_index()
    HLAs = open(file, 'r').readline().strip()
    ss = '\t' * len(rep_cols)
    HLAs = [i.replace('\t','') for i in HLAs.split(ss)]
    cols = [(i, "") for i in start_cols]
    for HLA in HLAs:
        for rep_col in rep_cols:
            cols.append((HLA, rep_col))
    cols.append(('Ave', '')); cols.append(('NB', ''))
    df.columns = pd.MultiIndex.from_tuples(cols)
    return(df)

def get_reverse_basepair(base_seq):
    return("".join([dic_pair_base[c] for c in base_seq[::-1]]))

def sort_str_num(List, pat = r'input/(\w+)-.*D(\d+)', pstr=1, pnum=2) -> list:
    return(sorted(List, key = lambda x: (re.search(pat, x).group(pstr), int(re.search(pat, x).group(pnum)))))

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def read_vcfINFO_to_df(vcffile, info_col=8, value_col=9):
    arr = []
    with open(vcffile, 'r') as f:
        for l in f:
            if l[0] == '#': continue
            l = l.strip().split('\t')
            flags = l[info_col].split(':'); values=l[value_col].split(':')
            values[6] = values[6].strip('%')
            dic = dict(zip(flags, values))
            arr.append(dic)
    df = pd.DataFrame(arr, dtype='float')
    return(df)

def trans_arr_col_row(in_arr):
    in_arr = np.array(in_arr); col_num = in_arr.shape[1]
    return([in_arr[:,i] for i in range(col_num)])

def Pearsonr(x, y):
    return scs.pearsonr(x, y)[0]

def triangle_square(lines):
    a, b, c = lines[0], lines[1], lines[2]
    p = (a + b + c) / 2
    return(np.sqrt(p*(p-a)*(p-b)*(p-c)))

skip_samples = []

def get_skip_samples(bad_pair_file):
    global skip_samples
    skip_samples = []
    tmp = []
    with open(bad_pair_file, 'r') as f:
        for l in f:
            skip_samples.append(l.strip())
            tmp.extend(l.strip().split('-VS-'))
        skip_samples.extend(tmp)

def skip_file(files):
    del_files = []
    for i in files:
        for s in skip_samples:
            if re.search(s, i):
                del_files.append(i)
                break
    for i in del_files:
        files.remove(i)

def get_coding_size(sample, qcpath):
    if re.search(r'-VS-', sample):
        sample = re.search(r'(.*)-VS-', sample).group(1)
        j_file = glob(f'{qcpath}/{sample}/*qualimm.json')[0]
        with open(j_file, 'r') as fo:
            jc = json.load(fo)
        cs = float(jc['Coding_size'].replace(',', '')) / 1000000
        return(cs)
    else:
        j_file = glob(f'{qcpath}/{sample}/*qualimm.json')[0]
        with open(j_file, 'r') as fo:
            jc = json.load(fo)
        cs = float(jc['Coding_size'].replace(',', '')) / 1000000
        return(cs)

def mdf2circos(infile, outfile):
    with open(infile, 'r') as f:
        region_num = defaultdict(dict)
        for l in f:
            l = l.strip().split('\t')
            start = (int(l[1])//1000000)*1000000
            end = start + 1000000
            region = f'{start}\t{end}'
            region_num[l[0]].setdefault(region, 0)
            region_num[l[0]][region] += int(l[2])
            
    with open(outfile, 'w') as fo:
        for Chr in region_num.keys():
            for k, v in region_num[Chr].items():
                fo.write(f'{Chr}\t{k}\t{v}\n')

def get_target_region(tar_region_file):
    with open(tar_region_file, 'r') as f:
        dic_region = defaultdict(list)
        for l in f:
            l = l.strip().split('\t')
            dic_region[l[0]].append([int(l[1]),int(l[2])])
    return(dic_region)

def mhc_in_region(mhc_site, tar_region_file):
    dic_region = get_target_region(tar_region_file)
    outfile = mhc_site + '.tarR_num'
    fo = open(outfile, 'w')
    with open(mhc_site, 'r') as f:
        dic_count = defaultdict(int)
        for l in f:
            l = l.strip().split('\t')
            p = int(l[1])
            for p1,p2 in dic_region[l[0]]:
                if p >= p1 and p <= p2:
                    dic_count[f'{l[0]}\t{p1}\t{p2}'] += int(l[2])
    
    for k,v in dic_count.items():
        fo.write(f'{k}\t{v}\n')
    fo.close()

def reB_pos(infile, outfile, Rnum):
    fo = open(outfile, 'w')
    with open(infile, 'r') as f:
        Chr = 'chrNone'
        start, end = 0, 0
        Rk = Rnum/1000000
        for l in f:
            l = l.strip().split('\t')
            dis = int(l[2]) - int(l[1])
            if Chr != l[0]:
                if Chr != 'chrNone':
                    print(f'chr - {Chr} {Chr} 0 {end} {Chr}')
                start = 1
                Chr = l[0]
            else:
                start = end + 1
            end = start + dis
            rpkm = (int(l[3])/Rk) / (dis/1000)
            #rpkm = int(l[3]) / (dis/1000)
            line = f'{Chr}\t{start}\t{end}\t{rpkm}\n'
            fo.write(line)
        else:
            print(f'chr - {Chr} {Chr} 0 {end} {Chr}')
    fo.close()

def get_target_region(tar_region_file):
    with open(tar_region_file, 'r') as f:
        dic_region = defaultdict(list)
        for l in f:
            l = l.strip().split('\t')
            dic_region[l[0]].append([int(l[1]),int(l[2])])
    return(dic_region)

def mhc_in_region(mhc_site, tar_region_file):
    dic_region = get_target_region(tar_region_file)
    outfile = mhc_site + '.tarR_num'
    fo = open(outfile, 'w')
    dic_count = {}
    for Chr in dic_region.keys():
        for p1,p2 in dic_region[Chr]:
            s = f'{Chr}\t{p1}\t{p2}'
            dic_count[s] = 0

    with open(mhc_site, 'r') as f:
        for l in f:
            l = l.strip().split('\t')
            p = int(l[1])
            for p1,p2 in dic_region[l[0]]:
                s = f'{l[0]}\t{p1}\t{p2}'
                if p >= p1 and p <= p2:
                    dic_count[s] += int(l[2])
    
    for k,v in dic_count.items():
        fo.write(f'{k}\t{v}\n')
    fo.close()
    
def read_refgene_txt(txt_file = '/mnt/nfs/database/hg19/ucsc/refgene.20171024.txt'):
    dic = {}
    gene_chr = {}
    with open(txt_file, 'r') as f:
        for l in f:
            if re.match('^#', l): continue
            l = l.split('\t')
            gname = l[-4]
            if dic.get(gname): continue
            starts = l[9].split(','); starts.remove(''); starts = map(int, starts)
            ends = l[10].split(','); ends.remove(''); ends = map(int, ends)
            dic[gname] = list(zip(starts, ends))
            gene_chr[gname] = l[2]
    return(dic, gene_chr)

def get_AAfile(AApath, pair):
    if glob(f'{AApath}/{pair}/*snv.sindel.AAchange.xls.t'):
        return(glob(f'{AApath}/{pair}/*snv.sindel.AAchange.xls.t')[0])
    elif glob(f'{AApath}/{pair}/*snv.sindel.AAchange.xls.bak'):
        return(glob(f'{AApath}/{pair}/*snv.sindel.AAchange.xls.bak')[0])
    elif glob(f'{AApath}/{pair}/*snv.sindel.AAchange_v0.xls'):
        return(glob(f'{AApath}/{pair}/*snv.sindel.AAchange_v0.xls')[0])
    else:
        return(glob(f'{AApath}/{pair}/*snv.sindel.AAchange.xls')[0])

def read_MutInfo_to_df(name, AA_file, vaf_file, hotspot=None):
    with open(vaf_file, 'r') as f:
        vaf = defaultdict(dict); f.readline()
        for l in f:
            l = l.strip().split('\t'); l[0] = l[0].replace('chr', '')
            key = f'{l[0]}\t{l[1]}'
            vaf[key]['VAF'], vaf[key]['Alt_Reads'] = float(l[2])*100, int(l[4])
            
    cols_order = ['Chr', 'Pos', 'Gene', f'm.DNA', f'c.Protein']
    cols_num = [0, 1, 5, 9, 10]        
    with open(AA_file, 'r') as f:
        arr = []; f.readline()
        for l in f:
            l = l.strip().split('\t')
            tmp = dict(zip(cols_order, map(lambda x: l[x], cols_num)))
            key = tmp['Chr'] + '\t' + tmp['Pos']
            tmp[f'VAF'], tmp[f'Alt_Reads'] = vaf[key]['VAF'], vaf[key]['Alt_Reads']
            if hotspot: tmp['hotspot'] = 'Y' if f"chr{tmp['Chr']}_{tmp['Pos']}" in hotspot else 'N'
            arr.append(tmp)
        cols_order.extend(['VAF', 'Alt_Reads'])
        if len(arr)==0: return(pd.DataFrame(columns=cols_order+['Sample','hotspot'] if hotspot else cols_order))
        df_out = pd.DataFrame(arr)[cols_order+['hotspot'] if hotspot else cols_order]
        df_out['Sample'] = name
    return(df_out)

def MHC_cal(mhcfile):
    with open(mhcfile, 'r') as f:
        dic_pos = {}
        num = 0
        for l in f:
            if re.match('^NeoRank|None', l): continue
            l = l.split('\t')
            pos = f'{l[9]}_{l[10]}'
            if dic_pos.get(pos):
                continue
            dic_pos[pos] = 1
            num += 1
    return(num)

def MHC_cal_clonal(mhcfile):
    with open(mhcfile, 'r') as f:
        dic_pos = {}
        num, clonal, subclonal = 0, 0, 0
        for l in f:
            if re.match('^NeoRank|None', l): continue
            l = l.strip().split('\t')
            pos = f'{l[9]}_{l[10]}'
            if dic_pos.get(pos):
                continue
            dic_pos[pos] = 1
            num += 1
            if l[-1] == 'clonal':
                clonal += 1
            elif l[-1] == 'subclonal':
                subclonal += 1
    return(num, clonal, subclonal)

def file_lineNum(file):
    with open(file) as f:
        for i, l in enumerate(f):
            pass
        return(i + 1)

def cancer_type_rec(c_name):
    c_dic = {'肺':'肺癌', '食管|食道':'食管癌','胃':'胃癌', '卵巢':'卵巢癌','宫颈':'宫颈癌', '结肠|直肠':'结直肠癌',
         '乳腺':'乳腺癌', '鼻咽':'鼻咽癌','膀胱':'膀胱癌', '肾':'肾癌','肝':'肝癌', '胰':'胰腺癌', '黑色素':'黑色素瘤',
            '胆管':'胆管癌', 'Lung':'肺癌', 'NSCLC':'肺癌', 'menanoma':'黑色素瘤', 'NPC':'鼻咽癌'}
    if len(re.findall('癌', c_name)) > 1:
        return('多种肿瘤')
    me = re.search(r'(\w+)癌', c_name)
    if me:
        c_in = me.group(1)
    else:
        c_in = c_name
    for i in c_dic.keys():
        if re.search(i, c_in):
            return(c_dic[i])
    else:
        return(c_name)

def get_MMR_list(list_file='/mnt/nfs/database/global/yucemed/report/MMR.list'):
    global MMRgene
    MMRgene = []
    with open(list_file, 'r') as f:
        for l in f:
            MMRgene.append(l.strip())

def if_have_MMR_mutation(AAchange_file):
    with os.popen(f'cut -f6 {AAchange_file}') as f:
        for l in f:
            if l.strip() in MMRgene:
                return(1)
        else:
            return(0)

def get_ascatngs(f_pattern):
    if glob(f_pattern):
        file = glob(f_pattern)[0]
        contt = os.popen(f'egrep "rho|psi" {file}').read()
        Purity = float(re.search(r'rho\s(\S+)', contt).group(1))
        Ploidy = float(re.search(r'psi\s(\S+)', contt).group(1))
        if Purity == 1:
            Purity, Ploidy = np.nan, np.nan
    else:
        Purity, Ploidy = np.nan, np.nan
    return(Purity, Ploidy)

def msi_value(f_Pattern):
    if glob(f_Pattern):
        msi_file = glob(f_Pattern)[0]
        msi = os.popen(f'sed "1d" {msi_file} |cut -f3').readline().strip()
        return(float(msi))
    else:
        return(np.nan)

def if_have_fusion(f_Pattern):
    if glob(f_Pattern):
        fusion_file = glob(f_Pattern)[0]
        if file_lineNum(fusion_file) > 1:
            return(1)
        else:
            return(0)
    else:
        return(np.nan)

def if_hlaloh(f_Pattern):
    if glob(f_Pattern):
        hla_file = glob(f_Pattern)[0]
        hla = open(hla_file, 'r').read()
        if re.search('阳性', hla):
            return('Positive')
        else:
            return('Negative')
    else:
        return(np.nan)

def if_msi(inval):
    if inval > 20:
        return('MSI')
    else:
        return('MSS')

def tnb_level(inval):
    if inval < tnb_low:
        return('TNB-Low')
    elif inval > tnb_high:
        return('TNB-High')
    else:
        return('TNB-Medium')

def tmb_level(inval):
    if inval < tmb_low:
        return('TMB-Low')
    elif inval > tmb_high:
        return('TMB-High')
    else:
        return('TMB-Medium')

def ch2en(ct):
    ch2en = {'肺癌':'Lung Cancer', '肺腺癌':'LUAD', '小细胞肺癌':'SCLC', '肺鳞癌':'SquaLC', '肝癌':'Liver Cancer', '食管癌':'Esophagus Cancer', '胃癌':'Stomach Cancer', 
             '结直肠癌':'Colorectal Cancer', '胰腺癌':'Pancreatic Cancer', '其他':'Others', '胆管癌':'Cholangiocarcinoma', '膀胱癌':'Bladder Cancer',
            '多种肿瘤':'Multiple tumors', '宫颈癌':'Cervical Carcinoma', '乳腺癌':'Mammary Cancer', '卵巢癌':'Oophoroma', '胆囊癌':'gallbladder carcinoma',
            '肾癌':'Renal Carcinoma', '黑色素瘤':'Melanoma', '鼻咽癌':'Nasopharyngeal Carcinoma', '软组织肉瘤':'soft tissue sarcoma'}
    if ct in ch2en:
        return(ch2en[ct])
    else:
        return('Others')

def get_ycOne_region(ycOne='/mnt/nfs/database/global/yucemed/panel/YConePlus_v1.3.call.bed'):
    ycone_region = defaultdict(list)
    with open(ycOne, 'r') as f:
        for l in f:
            l = l.strip().split('\t')
            ycone_region[l[0]].append([l[1], l[2]])
    return(ycone_region)

def get_panel_result(Chr, Pos):
    for p1,p2 in ycone_region[Chr]:
        if Pos >= p1 and Pos <= p2:
            return(1)
    else:
        return(0)

def TNB_exon_panel(mhcfile):
    with open(mhcfile, 'r') as f:
        dic_pos = {}
        tnb, panel = 0, 0
        for l in f:
            if re.match('^NeoRank|None', l): continue
            l = l.split('\t')
            pos = f'{l[9]}_{l[10]}'
            if dic_pos.get(pos):
                continue
            dic_pos[pos] = 1
            panel += get_panel_result(f'chr{l[9]}',l[10])
            tnb += 1
    return(tnb/30, panel/1.5)

def single_get_first(unicode1):
    str1 = unicode1.encode('gbk')
    try:
        ord(str1)
        return str1.decode('gbk')
    except:
        asc = str1[0] * 256 + str1[1] - 65536
        if asc >= -20319 and asc <= -20284:
            return 'A'
        if asc >= -20283 and asc <= -19776:
            return 'B'
        if asc >= -19775 and asc <= -19219:
            return 'C'
        if asc >= -19218 and asc <= -18711:
            return 'D'
        if asc >= -18710 and asc <= -18527:
            return 'E'
        if asc >= -18526 and asc <= -18240:
            return 'F'
        if asc >= -18239 and asc <= -17923:
            return 'G'
        if asc >= -17922 and asc <= -17418:
            return 'H'
        if asc >= -17417 and asc <= -16475:
            return 'J'
        if asc >= -16474 and asc <= -16213:
            return 'K'
        if asc >= -16212 and asc <= -15641:
            return 'L'
        if asc >= -15640 and asc <= -15166:
            return 'M'
        if asc >= -15165 and asc <= -14923:
            return 'N'
        if asc >= -14922 and asc <= -14915:
            return 'O'
        if asc >= -14914 and asc <= -14631:
            return 'P'
        if asc >= -14630 and asc <= -14150:
            return 'Q'
        if asc >= -14149 and asc <= -14091:
            return 'R'
        if asc >= -14090 and asc <= -13119:
            return 'S'
        if asc >= -13118 and asc <= -12839:
            return 'T'
        if asc >= -12838 and asc <= -12557:
            return 'W'
        if asc >= -12556 and asc <= -11848:
            return 'X'
        if asc >= -11847 and asc <= -11056:
            return 'Y'
        if asc >= -11055 and asc <= -10247:
            return 'Z'
        return ''

def getPinyin(string):
    if string == None:
        return None
    lst = list(string)
    charLst = []
    for l in lst:
        charLst.append(single_get_first(l))
    return ''.join(charLst)

def getGender(gender):
    if gender == 'M':
        return('XY')
    elif gender == 'F':
        return('XX')
    else:
        return('Unknown')
