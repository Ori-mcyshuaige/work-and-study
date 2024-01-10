# -*- coding: utf-8 -*-
from .InviteFriends import *
import requests, subprocess

##########################################################################
## RESTAPI functions
## get genome, transcript sequence, variation annotation if related info by ensembl RESTAPI

server = "http://grch37.rest.ensembl.org"


def RESTAPI_get_loc_with_transID(transID):
    """e.g transID='ENST00000330501' """
    ext = f"/overlap/id/{transID}?feature=transcript"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    for i in r.json():
        if i['transcript_id'] == transID:
            return(i['seq_region_name'], i['start'], i['end'], i['strand'])
    
def RESTAPI_get_seq_of_region(Chr, p1, p2, strand=1):
    """return a plain sequence text"""
    ext = f"/sequence/region/human/{Chr}:{p1}..{p2}:{strand}?"
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    return(r.text)

def RESTAPI_get_tanscriptSeq_with_transID(transID):
    """e.g transID='ENST00000330501', return a fasta text with region and base info"""
    ext = f"/sequence/id/{transID}?type=cds"
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    return(r.text)

def RESTAPI_get_vepanno_with_gvar(gvar):
    '''6:g.51889864_51890506del'''
    ext = f"/vep/human/hgvs/{gvar}?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    return(r.json())

def RESTAPI_get_vepanno_with_transvar(transvar):
    '''ENST00000366667:c.803C>T'''
    ext = f"/vep/human/hgvs/{transvar}?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    return(r.json())

def RESAPI_pubmedID_with_var(gvar):
    '''6:g.51889864_51890506del'''
    ext = f"/vep/human/hgvs/{gvar}?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        return([])
    
    for i1 in r.json():
        if 'colocated_variants' in i1:
            for i2 in i1['colocated_variants']:
                if 'pubmed' in i2:
                    return(i2['pubmed'])
    return([])

##########################################################################
## Synapse downlad
def get_from_synapse_id(synid, dest):
    import synapseclient
    syn = synapseclient.Synapse()
    syn.login()
    entity = syn.get(synid, downloadLocation=dest) # syn22153884
 
##########################################################################
## Run Command, execut shell command and return
def status_message(msg):
    print(msg)
    sys.stdout.flush()

def run_cmd(cmd, msg=''):
    status_message(cmd)
    if ',' in msg:
        begin, finish = msg.split(',')
        status_message(begin)
    else:
        finish = msg
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,)
    except subprocess.CalledProcessError as err:
        error_msg = 'Error happend!: {}\n{}'.format(err, err.output)
    else:
        error_msg = ''
    if not error_msg:
        status_message(finish)
        return True
    else:
        status_message(error_msg)
        return False

