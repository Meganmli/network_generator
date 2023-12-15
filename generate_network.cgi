#!/usr/local/bin/python3

import cgi, json, os, requests, sys
import pandas as pd  
from urllib.request import urlopen; from requests.adapters import HTTPAdapter, Retry

from intermine.webservice import Service    # intermine==1.13.0
import mygene                               # mygene==3.2.2
from biothings_client import get_client     # biothings-client==0.3.1
    
''' MouseGene class stores identifier information for each mouse gene '''  
class MouseGene():
    def __init__(self, mgiid, entrezid, symbol):
        # self mouse attributes
        self.mgiid = mgiid
        self.entrezid = entrezid
        self.symbol = symbol
        self.mgiid_interactants = list()
        self.input = "false"
        # orthologous human attributes
        self.human_entrezid = None
        self.human_symbol = None
        self.human_ortholog_score = None
        self.human_pubs = None
        self.human_disease_assoc = None
    
    # ----- Getters -----
    def get_mgiid(self): return self.mgiid  
    def get_entrezid(self): return self.entrezid  
    def get_symbol(self): return self.symbol  
    def get_human_entrezid(self): return self.human_entrezid  
    def get_human_symbol(self): return self.human_symbol  
    def get_human_ortholog_score(self): return self.human_ortholog_score  
    def get_human_pubs(self): return self.human_pubs  
    def get_input(self): return self.input
    def get_mgiid_interactants(self): return self.mgiid_interactants
    
    # ----- Setters -----
    ''' Sets input to true if this MouseGene is the input/seed gene '''
    def set_input_true(self): self.input = 'true'
    ''' Sets orthologous human entrezid '''
    def set_human_entrezid(self, entrezid): self.human_entrezid = entrezid
    ''' Sets orthologous human symbol '''
    def set_human_symbol(self, symbol): self.human_symbol = symbol
    ''' Sets orthologous human confidence score '''
    def set_human_ortholog_score(self, score): self.human_ortholog_score = score
    ''' Sets orthologous human gene publication count '''
    def set_human_publications(self, pubs): self.human_pubs = pubs
    
    # ----- Methods -----
    ''' Pretty prints this MouseGene object'''
    def describe_me(self):
        print(f"{self.mgiid} with entrezid {self.entrezid} and symbol {self.symbol} \
              \n\twith human ortholog (score {self.human_ortholog_score}): entrezid {self.human_entrezid} and symbol {self.human_symbol} \
              \n\twith {self.human_pubs} publications \
              \n\twith {self.mgiid_interactants} interactants")

    ''' Add interactant for this MouseGene '''
    def add_interactant(self, mgiid_interactant):
        self.mgiid_interactants.append(mgiid_interactant)
        
''' Given a gene_name, finds all protein-protein M. Musculus ineractions from AllianceMine''' 
def get_interactions(results, gene_name):
    mouse_genes_lst = list()
    
    # convert gene symbol name to MGI gene id
    service = Service('https://www.alliancegenome.org/alliancemine/')    
    query = service.new_query("Gene")
    query.add_constraint("symbol", "=", gene_name, code="A")
    query.add_constraint("organism", "LOOKUP", 'Mus musculus', code="B")

    mgi_gene_id = ''
    for row in query.rows(): 
        # print(row["primaryIdentifier"], row["name"], row["symbol"], row["secondaryIdentifier"])
        mgi_gene_id = row["primaryIdentifier"]

    # store species ids in list (ex: MGI:XX)
    df = pd.DataFrame(columns=['geneA_entrezid','geneA_symbol','geneB_entrezid','geneB_symbol'])
    # prepare url
    filter_interactorspecies = ('Mus musculus').replace(' ','%20')
    url = f'https://www.alliancegenome.org/api/gene/{mgi_gene_id}/interactions?filter.interactorSpecies={filter_interactorspecies}&filter.joinType=molecular_interaction'
    try: response = urlopen(url)
    except: return results, mouse_genes_lst
    data_json = json.loads(response.read())
    rows = data_json['results']

    geneA_id = ''; geneA_symbol = ''; geneA_entrezid = ''; seed_MouseGene = ''
    # loop through each geneB interactant
    for i, row in enumerate(rows):
        geneA_id = row['geneA']['id']; geneA_symbol = row['geneA']['symbol']
        geneA_entrezid, geneA_symbol  = convert_to_entrez_symbol(geneA_id, geneA_symbol)
        
        # only attempt finding geneB if geneA has entrezid
        if geneA_entrezid != "None":
            # add gene A once
            if i == 0:
                results['interactants'].append({'mouse_id': geneA_id, 'mouse_entrezid': geneA_entrezid, 'input': 'true',
                                                'human_entrezid':'','human_symbol':'','score':'0','pubs':'0'}) 
                seed_MouseGene = MouseGene(geneA_id, geneA_entrezid, geneA_symbol)
                seed_MouseGene.set_input_true()
                
            geneB_id = row['geneB']['id']; geneB_symbol = row['geneB']['symbol']; geneB_entrezid = ''
            geneB_entrezid, geneB_symbol  = convert_to_entrez_symbol(geneB_id, geneB_symbol)
            
             # only add if geneB has entrezid
            if geneB_entrezid != "None":
                mouse_genes_lst.append(MouseGene(geneB_id, geneB_entrezid, geneB_symbol))
                seed_MouseGene.add_interactant(geneB_id)
    
                # add gene Bs
                results['interactants'].append({'mouse_id': geneB_id, 'mouse_entrezid': geneB_entrezid,'input': 'false',
                                                'human_entrezid':'','human_symbol':'','score':'0','pubs':'0'})
                # add gene A to gene B interaction
                results['interactions'].append({'geneA_interactant': geneA_id, 'geneB_interactant': geneB_id})
            
    # only append if MouseGene initialized
    if seed_MouseGene != '': mouse_genes_lst.append(seed_MouseGene)
    return results, mouse_genes_lst

''' Converts a given query symbol to entrez ID, if present in bioclient or mygene'''
def convert_to_entrez_symbol(query, query_symbol):
    '''
    INPUT:
        query: str of gene database id
        query_symbol: str of gene symbol 
    OUTPUT:
        str of ({entrez gene id}, {gene symbol})
    '''
    # set mouse scope and taxon id
    scope = "MGI"; species_id = "10090"; fields = ['query,symbol,entrezgene']
    mg = get_client('gene')
    
    # (1) try bioclient
    conversions = mg.getgene(query, fields=fields, species=species_id, verbose=False)   
    # (1a) successful 
    if isinstance(conversions, list):
        conversions = conversions[-1] # return last item in list, most recent version
        return conversions['entrezgene'], conversions['symbol']
    # (1b) successful 
    elif (isinstance(conversions, dict)) and ('entrezgene' in conversions):
        return conversions['entrezgene'], conversions['symbol']
    
    # (2) fails come here, bioclient again using symbol
    conversions = mg.querymany([query_symbol], scopes='symbol', fields=fields, species=species_id, verbose=False)
    # (2a) successful 
    if isinstance(conversions, list):
        conversions = conversions[-1] # return last item in list, most recent version
        if 'entrezgene' in conversions:
            return conversions['entrezgene'], conversions['symbol']
                        
    # (3) fails come here, try mygene method
    mg = mygene.MyGeneInfo()
    try:
        conversions = mg.querymany([query], scopes=scope, fields=fields, species=species_id, verbose=False)[0]
    except requests.exceptions.ConnectionError:
        pass
    # (3a) failed
    if 'notfound' in conversions:
        # (4) mygene method using symbol
        conversions = mg.querymany([query_symbol], scopes='symbol', fields=fields, species=species_id, verbose=False)   
        if 'notfound' in conversions:
            return 'None','None'
        
    # (3b or 4b) successful        
    try: return conversions['entrezgene'], conversions['symbol']
    except: return 'None','None'
   
''' Handler that collects ortholog for each input gene '''
def get_orthologs(results, MouseGenes):
    '''
    INPUT:
       results: dict
       MouseGenes: lst of MouseGene objects
    OUTPUT: 
       results: dict with orthologs
       MouseGenes: lst of MouseGene objects with orthologs
   '''
    # for each entrezid. find best human ortholog
    for MouseGene in MouseGenes:
        mouse_entrezid = MouseGene.get_entrezid()
        entrezid, symbol, score = find_best_ortholog_diopt(mouse_entrezid)
        MouseGene.set_human_entrezid(entrezid)
        MouseGene.set_human_symbol(symbol)
        MouseGene.set_human_ortholog_score(score)

    for interactant_dict in results['interactants']:
        mouse_entrezid = interactant_dict['mouse_entrezid']
        entrezid, symbol, score = find_best_ortholog_diopt(mouse_entrezid)
        interactant_dict['human_entrezid'] = entrezid
        interactant_dict['human_symbol'] = symbol
        interactant_dict['score'] = str(score)
    return results, MouseGenes
    
''' From DIOPT, finds best ortholog(s) with highest diopt score and highest diopt confidence '''
def find_best_ortholog_diopt(entrez_gene_id: str):
    assert isinstance(entrez_gene_id, str)
    
    # setup session with 
    s = requests.Session()
    retries = Retry(total=10,
                    # backoff_factor 0.1, sleep() will sleep for [0.05s, 0.1s, 0.2s, 0.4s, ...] between retries
                    backoff_factor=0.1,
                    # force a retry if the status code returned is 500, 502, 503 or 504
                    status_forcelist=[ 500, 502, 503, 504 ])
    s.mount('http://', HTTPAdapter(max_retries=retries))
    endpoint = f'https://www.flyrnai.org/tools/diopt/web/diopt_api/v9/get_orthologs_from_entrez/10090/{entrez_gene_id}/9606/exclude_score_less_2'
    r = s.get(endpoint, headers={"Content-Type" : "application/json"})
    
    if not r.ok:
      r.raise_for_status()
      sys.exit()
 
    decoded = r.json(); search_results = []
    search_results.append(decoded)
   
    clean_results_dict = {}
    for search_result in search_results:
      # Extract results dictionary 
      results_dict = search_result['results']
      search_details_dict = search_result['search_details']
      if not results_dict: continue 
      else:
        # Extract input Entrez ID's results dictionary
        entrez_gene_id = list(results_dict.keys())[0]   
        # Collect model organism results
        for match_species in results_dict[entrez_gene_id].keys():
          # Add input search/query info
          search_details_dict_subset = {}
          search_details_dict_subset['input_species'] = search_details_dict['input_species']
          search_details_dict_subset['input_entrez_geneid'] = search_details_dict['gene_details'][0]['geneid']
          search_details_dict_subset['input_symbol'] = search_details_dict['gene_details'][0]['symbol']
          # Use update method to move search details to front
          match_species_dict = results_dict[entrez_gene_id][match_species]
          search_details_dict_subset.update(match_species_dict)
          clean_results_dict[match_species] = search_details_dict_subset
          
    # no overall results, empty
    if not clean_results_dict: return "none", "none", 0
    
    diopt_confidence = {"high":3, "moderate":2, "low":1}
    # pointers to track current best scores
    best_score = 0; best_confidence_score = 0; best_score_count = 0
    # final best ortholog
    human_entrezid = 'none'; human_symbol = 'none'; final_confidence_score = 0
    for k, v in clean_results_dict.items():
        t_best_score = v["score"]
        t_confidence_score = diopt_confidence[v["confidence"]]
        t_best_score_count = v["best_score_count"]
        # first tie breaker
        if t_best_score > best_score:
            best_score = t_best_score; best_confidence_score = t_confidence_score; best_score_count = t_best_score_count
            human_entrezid = v["geneid"]; human_symbol = v["symbol"]; final_confidence_score = diopt_confidence[v["confidence"]]
        # second tie breaker
        elif t_best_score == best_score:
            if t_confidence_score > best_confidence_score:
                best_confidence_score = t_confidence_score; best_score_count = t_best_score_count
                human_entrezid = v["geneid"]; human_symbol = v["symbol"]; final_confidence_score = diopt_confidence[v["confidence"]]
            elif t_confidence_score == best_confidence_score:
                if t_best_score_count > best_score_count:
                    best_score_count = t_best_score_count
                    human_entrezid = v["geneid"]; human_symbol = v["symbol"]; final_confidence_score = diopt_confidence[v["confidence"]]
    return human_entrezid, human_symbol, final_confidence_score
        
''' Retrieves publication counts from entrezid '''
def get_publications(results, MouseGenes):
   '''
   INPUT:
       results: dict
       MouseGenes: lst of MouseGene objects
   OUTPUT: 
       results: dict with publication counts
       MouseGenes: lst of MouseGene objects with publication counts
   '''
   generefs = pd.read_csv("./data/publications.csv")
   entrezid_to_refs_dict = dict(zip(generefs['entrezid'], generefs['PubMed_refs']))

   for MouseGene in MouseGenes:
        human_entrezid = MouseGene.get_human_entrezid()
        pubs = entrezid_to_refs_dict[human_entrezid]
        if (pubs == 0) or (pubs != pubs): pubs = "1"
        MouseGene.set_human_publications(str(pubs))
   
   for interactant_dict in results['interactants']:
       human_entrezid = interactant_dict['human_entrezid']
       pubs = entrezid_to_refs_dict[human_entrezid]
       if (pubs == 0) or (pubs != pubs): pubs = "1"
       interactant_dict['pubs'] = str(pubs)
   return results, MouseGenes

def main():
    # DEBUGGING EDGE CASES
    # print(find_best_ortholog_diopt('-1'))
    # sys.exit()
    
    # site launched, ready to go!
    print("Content-Type: application/json\n\n")
    
    # retrieve user input
    form = cgi.FieldStorage()
    term = form.getvalue('search_term')
    
    # term = "Rap1gds1" # "Loricrin" 
     
    results = { 'match_count': "0", 'interactants': list(), 'interactions': list() }
    # (1) get protein-protein interactions
    results, MouseGenes = get_interactions(results, term)
    # (1a) catch when no interactions found (ex: Glyatl3)
    if len(MouseGenes) == 0: 
        results = {"message":"No interactions found for this gene.\nPlease input another gene!"}
        print(json.dumps(results))
        
    # (2) get orthologs and confidence scores
    results, MouseGenes = get_orthologs(results, MouseGenes)
    
    # (3) get publication counts
    results, MouseGenes = get_publications(results, MouseGenes) 
    
    # (4) prepare results dictionary in json format using MouseGene object to build results last
    results = { "match_count": "0", "interactants": list(), "interactions": list() }
    for M in MouseGenes:
        # M.describe_me()
        mouse_id = str(M.get_mgiid())
        input = M.get_input()
        human_entrezid = str(M.get_human_entrezid())
        human_symbol = M.get_human_symbol()
        pubs = str(M.get_human_pubs())
        score = str(M.get_human_ortholog_score())
        
        results["interactants"].append({"mouse_id":mouse_id, "input":input,
                                        "human_entrezid":human_entrezid,"human_symbol":human_symbol,
                                        "pubs":pubs, "score":score})
        if input == 'true':
            geneB_mouseids = M.get_mgiid_interactants()
            for geneB_mouseid in geneB_mouseids:
                results["interactions"].append({"geneA_interactant":mouse_id, "geneB_interactant":geneB_mouseid})
            results["match_count"] = len(geneB_mouseids)
      
    # dump results to search.js
    print(json.dumps(results))

if __name__ == '__main__':
    main()
