from HappyWorkFriends.InviteFriends import *
from pymongo import MongoClient

def get_sample_info(s, items):
    tmp = {}
    for i in items:
        try:
            if i == 'QCbed':
                tmp[i] = os.path.basename(s[i])
            else:
                tmp[i] = s[i]
        except KeyError:
            tmp[i] = np.nan
    else:
        return(tmp)

def get_somatic_info(s, AA=False):
    if not s.get('DNAsomatic'): return([],np.nan)
    tmp = []; AA = 0
    for k in s['DNAsomatic']:
        if AA:
            if not k['Impact'] in ['HIGH','MODERATE','Impact']: continue
        sm = get_sample_info(k,['Chr','Pos','Gene','Exon','Ref','Alt','Type','Impact'])
        if not sm: continue
        sm['VAF'] = k.get('VAF') or np.nan
        sm['bed'] = os.path.basename(s['QCbed'])
        try:
            sm['SampleID'],sm['Contact'],sm['Ctype'] = s['SampleID'],s['Contact'],s['Ctype']
        except KeyError:
            pass
        tmp.append(sm); AA += 1
    TMB = AA/(int(s['Coding'])/1000000)
    return(tmp, TMB)

def db_query(db, cond={"SampleID":re.compile('XJ')}, dpth=500, Somatic=False, NeoAntigen=False, Fusion=False, Rela=False,
             items=['SampleID','Contact','QCbed','Ctype','Depth','Workdir']):
    """if Rela=True, Somatic=True is required"""
    if Rela == True: Somatic=True
    rslt = []; somatic = []; NeoAnti = []; fusion = []; rela = []; rel_sm = []
    for s in db.find(cond):
        if s['PatientF'] == 1:
            tmp = get_sample_info(s, items)
            if not tmp: continue
            if tmp['Depth'] < dpth: continue
            AA = 0; relaAA = 0
            if Somatic:
                sm, AA = get_somatic_info(s)
                if s['Coding'] == 0: print(s['SampleID'])
                tmp['TMB'] = AA
                somatic.extend(sm)

            if Fusion:
                if s.get('DNAfusion'):
                    for k in s['DNAfusion']:
                        fusion.append(k)

            if Rela:
                me=re.search(r'^(\w{8})',tmp['SampleID'])
                me2=re.search(r'^(\w{10})',tmp['SampleID'])
                if me:
                    f = me.group(1); f2 = me2.group(1)
                    for s in db.find({"SampleID":re.compile(f)}):
                        if not tmp['SampleID'] in s['SampleID'] and not f2 in s['SampleID']:
                            d = {}; sm, relaAA = get_somatic_info(s)
                            if not sm: continue
                            d['sample1'],d['sample2'],d['Ctype'] = tmp['SampleID'],s['SampleID'],tmp['Ctype']
                            d['dir1'],d['dir2'] = tmp['Workdir'],s['Workdir']
                            d['tmb1'],d['tmb2'] = AA, relaAA
                            rela.append(d); rel_sm.extend(sm)
            rslt.append(tmp)
    return(pd.DataFrame(rslt).dropna(subset=['SampleID']).drop_duplicates(subset='SampleID', keep='last'),
           pd.DataFrame(somatic).drop_duplicates(subset=['SampleID','Chr','Pos'], keep='last'),
           NeoAnti,fusion,
           pd.DataFrame(rela).drop_duplicates(subset=['sample1','sample2'], keep='last'),
           pd.DataFrame(rel_sm).dropna(subset=['SampleID']).drop_duplicates(subset=['SampleID','Chr','Pos'], keep='last'))
