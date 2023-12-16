#!/usr/local/bin/python3
import cgi, json
import os

import psycopg2 # psycopg2-binary==2.9.9

def main():
    print("Content-Type: application/json\n\n")
    # retrieve user input
    form = cgi.FieldStorage()
    # jQuery sends user-typing input as "term"
    term = form.getvalue('term') 

    # update connection credentials
    username = 'UPDATE ME'; password = 'UPDATE ME'
    host = 'mgi-adhoc.jax.org'; port = 5432; db_name = 'mgd'; jdbc_url = "jdbc:postgresql://mgi-adhoc.jax.org:5432/mgd "
    conn = psycopg2.connect(database=db_name, host=host, user=username, password=password, port=port)
    cursor = conn.cursor()
    
    # execute query
    # query 10 gene names from [mrk_label] where their (_marker_key) maps to (term) 'protein coding gene' in [mrk_mcv_cache]
    qry = """
        SELECT DISTINCT label as name FROM mrk_label 
        WHERE LOWER(label) LIKE LOWER(%s)
        AND labeltype LIKE 'MS'
        AND labeltypename LIKE 'current symbol'
        AND _marker_key IN 
            (SELECT _marker_key FROM mrk_mcv_cache 
            WHERE term LIKE 'protein coding gene')
        LIMIT 10
    """ 
    cursor.execute(qry, ( '%' + term + '%', ))
    
    # organize results and count matches
    results = []
    for gene_name in cursor:
        results.append({'value':gene_name, 'label':gene_name})

    # close connection
    conn.close()
    
    # dump results to search.html
    print(json.dumps(results))

if __name__ == '__main__':
    main()
