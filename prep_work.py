#!/usr/local/bin/python3

import os 
import pandas as pd  

''' RUN ONCE. Reads in gene2pubmed ftp dump and collapses it into unique genes '''
def create_publications_file():
    # Source doc: https://nam02.safelinks.protection.outlook.com/?url=https%3A%2F%2Fftp.ncbi.nlm.nih.gov%2Fgene%2FDATA%2FREADME&data=05%7C01%7Cmli186%40jhu.edu%7C774b48caac7d4f472f7408db9855ef94%7C9fa4f438b1e6473b803f86f8aedf0dec%7C0%7C0%7C638271267460540386%7CUnknown%7CTWFpbGZsb3d8eyJWIjoiMC4wLjAwMDAiLCJQIjoiV2luMzIiLCJBTiI6Ik1haWwiLCJXVCI6Mn0%3D%7C3000%7C%7C%7C&sdata=0zDloSt9Sn1gjzsm%2Btk6uasLV45YF%2BYoNGxp7eSho1U%3D&reserved=0
    url='https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz' 
    gene2pubmed = pd.read_csv(url, sep='\t') 
    gene2pubmed = gene2pubmed.sort_values(by='GeneID')
    gene2pubmed = gene2pubmed.rename(columns={'GeneID': 'input_entrez_geneid'})
    
    # Store number of references per gene in df
    generefs = gene2pubmed.groupby("input_entrez_geneid")["PubMed_ID"].count().reset_index(name="PubMed_refs")
    entrezid_to_refs_dict = dict(zip(generefs['input_entrez_geneid'], generefs['PubMed_refs']))
    
    df = pd.DataFrame(entrezid_to_refs_dict.items(), columns=['entrezid', 'PubMed_refs'])
    df.to_csv("./data/publications.csv", index=False) 
    
    
def main():
    # BEFORE site launch, generate publications file 
    create_publications_file()
    

if __name__ == '__main__':
    main()
