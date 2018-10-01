#! /usr/bin/python


import shlex, subprocess
import sys


# TODO 
# try except if path, file, format
# 


class Run:
    def __init__(self, seqdata, biodb, add = None):

        # takes ID of genom goes in to bioSQL and get it runs 
        from BioSQL import BioSeqDatabase
        username = raw_input("Please enter user name: ")
        password = raw_input("and password: ")

        server = BioSeqDatabase.open_database(driver = "psycopg2", user = username, passwd = password, host = "dbpg-ifi-utv.uio.no", db = "rnammer")

        db = server[biodb]             # should be possible to input


        #seq_record = db.lookup(gi= seqdata).format("fasta") # extraxt string of only alphabet acgttgca
        
        name    = seqdata
        acc_id  = name.split(".")[0]
        version = name.split(".")[1]

        print acc_id 
        
        seq_record = db.lookup(accession= acc_id).format("fasta") # extraxt string of only alphabet acgttgca
        


        fasta_file_ish = seq_record # "%s \n%s" %(first_line,seq_record)

        f = open("tmpRnammmerFasta.txt","w")
        f.write(fasta_file_ish)
        f.close()
        
        # TODO check if genbank ID or fasta file
        args_str = "perl rnammer/rnammer -S bac -m lsu,ssu,tsu -gff - "
        
        args = shlex.split(args_str + "tmpRnammmerFasta.txt") 
        
        with open('r2.txt', 'wb', 0) as file:
            subprocess.call(args, stdout=file)

x = Run("NC_006462.1", "refseq")
print x

