from .InviteFriends import *
from .Olivia import get_coding_size
from pymongo import MongoClient

def read_vardict_muts(sample, mut_json):
    dic_rename = {'chromosome':'Chr', 'start':'Pos',
                  'variant_freq':f'VAF.{sample}', 'variant_reads':f'Alt_Reads.{sample}'}
    with open(mut_json, 'r') as f:
        mut_js = json.load(f)
    df_mut = pd.DataFrame(mut_js['somatic_snvindel']); df_mut = df_mut[df_mut.impact.isin(['MODERATE','HIGH'])]
    df_mut.start = df_mut.start.astype(str); df_mut.variant_freq = df_mut.variant_freq * 100
    # strip p. symbol
    df_mut['c.Mut'] = [i[2:] for i in df_mut.hgvsp_abbr]
    
    cols = ['chromosome', 'start', 'variant_freq', 'variant_reads', 'c.Mut']
    df_mut = df_mut[cols].rename(columns=dic_rename)
    #del pos specify
    df_mut.loc[df_mut.Pos == '55242465', 'Pos'] = '55242464'
    df_mut.loc[df_mut.Pos == '55242466', 'Pos'] = '55242465'
    return(df_mut)


def read_MutInfo_to_df(name, AA_file, vaf_file, hotspot=None):
    with open(vaf_file, 'r') as f:
        vaf = defaultdict(dict); f.readline()
        for l in f:
            l = l.strip().split('\t')
            key = f'{l[0]}:{l[1]}'
            vaf[key]['VAF'], vaf[key]['Alt_Reads'] = float(l[2])*100, int(l[4])
            
    cols_order = ['Chr', 'Pos', 'Gene', f'm.DNA', f'c.Protein']
    cols_num = [0, 1, 5, 9, 10]        
    with open(AA_file, 'r') as f:
        arr = []; f.readline()
        for l in f:
            l = l.strip().split('\t')
            tmp = dict(zip(cols_order, map(lambda x: l[x], cols_num)))
            tmp['Chr'] = 'chr' + tmp['Chr']
            key = tmp['Chr'] + ':' + str(tmp['Pos'])
            tmp[f'VAF'], tmp[f'Alt_Reads'] = vaf[key]['VAF'], vaf[key]['Alt_Reads']
            if hotspot: tmp['hotspot'] = 'Y' if f"chr{tmp['Chr']}_{tmp['Pos']}" in hotspot else 'N'
            arr.append(tmp)
        cols_order.extend(['VAF', 'Alt_Reads'])
        if len(arr)==0: return(pd.DataFrame(columns=cols_order+['Sample','hotspot'] if hotspot else cols_order))
        df_out = pd.DataFrame(arr)[cols_order+['hotspot'] if hotspot else cols_order]
        df_out['Sample'] = name
    return(df_out)

def missing_var_info(Chr, Pos, vcf_files, strict_files, info_col=8, value_col=9):
    """find out the missing variant information"""
    tmp = {}
    for vcf_file in vcf_files:
        with open(vcf_file, 'r') as f:
            for l in f:
                l = l.split('\t')
                if l[0] == Chr and l[1] == Pos:
                    tmp['Chr'], tmp['Pos'] = Chr.replace('chr',''), Pos
                    flags = l[info_col].split(':'); values=l[value_col].split(':')
                    tmp['VAF'] = values[6].strip('%')
                    tmp['Alt_Reads'] = int(values[5])
                    break
    for strict_file in strict_files:
        with open(strict_file, 'r') as f:
            for l in f:
                l = l.split('\t')
                if l[0] == Chr and l[1] == Pos:
                    tmp['VAF'] = tmp['VAF'] + f'({l[4].strip()})'
                    break
    return(tmp)

def add_sample_merged_column(writer, sheet, samples, start_col):
    alphabets = string.ascii_uppercase; alphabets = list(alphabets) + [f'A{i}' for i in alphabets]
    workbook = writer.book
    worksheet = writer.sheets[sheet]
    merge_format = workbook.add_format({
    'bold': 1,
    'border': 1,
    'align': 'center'})
    
    samples = [re.match('(.*)\..*', i).group(1) for i in samples if '.' in i]
    samples = np.unique(samples)
    for s in samples:
        Rng = f'{alphabets[start_col]}1:{alphabets[start_col+2]}1'
        worksheet.merge_range(Rng, s, merge_format)
        start_col += 3

def get_standard_info(kind='GW001'):
    """kind: 'HD728','GW','cocktail_III','cocktail_VI' default:HD"""
    std_file_path = '/home/danxu/2.files/Test/Files/StandardReference/'
    std_file = std_file_path + kind + '.muts'
    df_std = pd.read_csv(std_file, sep='\t', names = ['Chr','Pos','Ref','Alt','Gene','c.Mut','ExpectP','type'], dtype='str')
    df_std = df_std[df_std.type == 'Somatic']; df_std.drop('type', axis=1, inplace=True)
    df_cmp = df_std[True ^ df_std['c.Mut'].isin(['Fusion','CNV'])].copy()
    fusion_genes = [i.split(':')[1] for i in df_std[df_std['c.Mut'] == 'Fusion'].Gene]
    cnv_genes = [re.search('(\w+)-', i).group(1) for i in df_std[df_std['c.Mut'] == 'CNV'].Gene]
    
    return(df_cmp, fusion_genes, cnv_genes)

def ordered_cols(in_cols):
    cols = ['Raw_Base(G)','Clean_Base(G)','Insert_size','Duplication_rate','Capture_rate','Depth_in_target','Depth/Raw_Base(G)',
            'Base_mapping_in_target/Raw_base(%)','Target_coverage',
       'Target_1000X','Target_500X','Target_300X','Target_100X','Target_10X','Clean_base_ratio','Mapping_rate','Mapping_quality',
       'Q20','Q30','GC','GC-AT_Seperation','N_Rate','Average_read_length','Read_length_stddev','Average_base_quality','Coding_size']
    in_cols = list(set(in_cols) & set(cols))
    d = {}; n=0
    for i in cols: d[i] = n; n+=1
    return(sorted(in_cols, key = lambda x: d[x]))

def conn_mongodb(func):
    def wrapper(*args, **kwargs):
        global conn
        conn = MongoClient('192.168.11.101', 27017)
        db = conn['admin']
        db.authenticate('medstat', 'MED0325$', mechanism = 'SCRAM-SHA-1')
        result = func(*args, **kwargs)
        conn.close()
        return(result)
    return(wrapper)

def conv2float(inval):
    inval = inval.strip('% Xbp'); inval = float(inval)
    return(inval)

def _read_qc_info(sample, qcfiles, qc_info):
    if sample in qc_info:
        sample = sample + '_mirror'
    for qcfile in qcfiles:
        with open(qcfile, 'r') as jf:
            qc_js = json.load(jf)
        for k,v in qc_js.items():
            if k == 'Coding_size':
                qc_info[sample][k] = float(v.replace(',', ''))
            else:
                qc_info[sample][k] = conv2float(v[1])
    if len(qcfiles) == 2:
        qc_info[sample]['Clean_Base(G)'] = round(qc_info[sample]['Clean_Base']/1000000000, 3)
        qc_info[sample]['Raw_Base(G)'] = round(qc_info[sample]['Clean_Base(G)']/(qc_info[sample]['Clean_base_ratio']/100), 3)
        qc_info[sample]['Depth/Raw_Base(G)'] = round(qc_info[sample]['Depth_in_target']/qc_info[sample]['Raw_Base(G)'], 2)
        qc_info[sample]['Base_mapping_in_target/Raw_base(%)'] = \
        round(qc_info[sample]['Depth_in_target']*qc_info[sample]['Coding_size']/(qc_info[sample]['Clean_Base']/(qc_info[sample]['Clean_base_ratio']/100))*100, 2)

@conn_mongodb
def read_qc_info(sample, path='/mnt/cfs/med18b/med_rd/danxu', opath='/mnt/cfs/med18b/med_rd/danxu/History_Stat', history=0):
    qc_info = defaultdict(dict)
    if glob(f'{path}/qc/{sample}/{sample}*json'):
        qcfiles = glob(f'{path}/qc/{sample}/{sample}*json')
        _read_qc_info(sample, qcfiles, qc_info)
    else:
        print(f'{sample} not exists in "{path}"')
    
    if history == 1:
        #me = re.search(r'^(\w{10})', sample)
        me = re.search(r'^DN\w(9\w{7})', sample)
        if me:
            flag = me.group(1)
            if glob(f'{opath}/*/qc/{flag}*'):
                qcfiles = glob(f'{opath}/*/qc/{flag}*/*.json')
                osample =  os.path.basename(os.path.dirname(fqcheck_file))
                _read_qc_info(osample, qcfiles, qc_info)
            else:
                db = conn.clindata
                my_set = db.Sinfo
                sinfos = my_set.find({"SampleID":{'$regex':flag}}).sort('_id',-1).limit(1)
                for sinfo in sinfos:
                    opath = sinfo['Workdir']
                    osample = sinfo['SampleID']
                    qcfiles = glob(f'{opath}/qc/{osample}/{osample}*json')
                    _read_qc_info('db.' + osample, qcfiles, qc_info)
        
    return(qc_info)

@conn_mongodb
def read_DBresult(sample):
    my_set = conn.clindata.Sinfo
    sinfo = my_set.find({"SampleID":sample})
    rslt = defaultdict(dict)
    if sinfo:
        N = 0
        for i in sinfo:
            if N == 0:
                flag = sample + '_manually'
                rslt[flag]['TMB'] = i['TMB']
                rslt[flag]['TNB'] = i['TNB'] if i.get('TNB') else np.nan
                rslt[flag]['MSI'] = i['MSI'] if i.get('MSI') else np.nan
            elif N == 1:
                rslt[sample]['TMB'] = i['TMB']
                rslt[sample]['TNB'] = i['TNB'] if i.get('TNB') else np.nan
                rslt[sample]['MSI'] = i['MSI'] if i.get('MSI') else np.nan
            else:
                break
            N += 1
        if N == 2:
            return(rslt)

def gen_fq_percent(fqfiles, odir, prefix):
    m_arr = np.linspace(0.5, 0.9, 5)[::-1]
    for p in m_arr:
        prefix = prefix + '_' + str(p)
        odir2 = f'{odir}/{prefix}'
        os.makedirs(f'{odir2}', exist_ok=True)
        b_name1 = os.path.basename(fqfiles[0]); b_name2 = os.path.basename(fqfiles[1])
        ofq1 = f'{odir2}/{b_name1}'; ofq2 = f'{odir2}/{b_name2}'
        fo1 = gzip.open(ofq1, 'wb'); fo2 = gzip.open(ofq2, 'wb')
        with gzip.open(fqfiles[0], 'r') as f1:
            with gzip.open(fqfiles[1], 'r') as f2:
                while True:
                    try:
                        lines1 = [f1.readline() for i in range(4)]
                        lines2 = [f2.readline() for i in range(4)]
                        if np.random.rand() <= p:
                            for l in lines1: fo1.write(l)
                            for l in lines2: fo2.write(l)
                    except EOFError:
                        break
        fo1.close(); fo2.close()

def get_hotspot_list(f_hotspot = '/mnt/nfs/database/global/yucemed/snv/hg19.hotspot.S30.bed'):
    hotspot = []
    with open(f_hotspot, 'r') as f:
        for l in f:
            l = l.split('\t')
            site = f'chr{l[0]}_{l[2]}'
            hotspot.append(site)
    return(hotspot)
