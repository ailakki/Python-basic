import fnmatch
import os
from Bio import SeqIO
from Bio import GenBank
from Bio.GenBank import utils

import re

def loadDB(catalog):
    from BioSQL import BioSeqDatabase
    import sys
    
    username = raw_input("Please enter user name: ")
    password = raw_input("and password: ")

    host = "dbpg-ifi-utv.uio.no"
    db_name = "rnammer"

    server = BioSeqDatabase.open_database(driver="psycopg2", user=username,passwd=password, 
            host=host, db=db_name)
    
    biodb_name = "empty"     # genebank problem ? se staving

    db  = "nodb"

    gi_rep = 1
    
    for gbff in catalog:
                               #server.remove_database(source)
        print gi_rep
        print gbff


        parser = GenBank.FeatureParser()
        #record = parser.parse(open(gbff))
        #records = SeqIO.parse(open(gbff),'genbank')
        records = GenBank.Iterator(open(gbff), parser)
        
        for x in records:
            if re.search("plasmid",x.description, re.IGNORECASE):
                continue
            print "Record name:"
            print x.id
            #print dir(x)

            if "Proteobacteria" == x.annotations["taxonomy"][1]:
                print x.annotations["taxonomy"][1]
                print x.annotations["taxonomy"][2]
                biodb_name = x.annotations["taxonomy"][2]
            else :
                print x.annotations["taxonomy"][1]
                biodb_name = x.annotations["taxonomy"][1]
            while True : 
                try :
                    db = server[biodb_name] 
                    #print "here"
                    break
                except KeyError :
                    #print ("Cannot find biodatabase with name %r making it" % source)
                    server.new_database(biodb_name)
                    server.commit()
            db.load([x])
        #record.annotations["gi"] = gi_rep 
        #print type(records)

        #print record.id
        gi_rep = gi_rep + 1

        #db.load([records])

    server.adaptor.commit()


matches = []
for root, dirnames, filenames in os.walk('/uio/hume/student-u34/ailaka/master/refseq/bacteria'):
    for filename in fnmatch.filter(filenames, '*.gbff'):
        matches.append(os.path.join(root, filename))
#print matches


loadDB(matches)

#loadDB(theSix,"theSixDB")
