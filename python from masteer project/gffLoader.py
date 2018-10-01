from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC

from BioSQL import BioSeqDatabase
import sys

username = raw_input("Please enter user name: ")
password = raw_input("and password: ")

server = BioSeqDatabase.open_database(driver = "psycopg2", user = username, passwd = password, host = "dbpg-ifi-utv.uio.no", db = "rnammer")


def gffLoader(seqdata,gffFile):
    
    # gloabal default
    biodb_name   = "biodb_name"  # dummy name
    db = "no_biodb"

    # Step 1 making secRecord

    name    = seqdata
    acc_id  = name.split(".")[0]
    version = name.split(".")[1]
    description = "Genome description and info" 

    #gi_id = ??

    S= Seq("MAAGVKQLADDRTLLDG", IUPAC.protein)

    rec = SeqRecord(
            S, 
            id= acc_id, 
            name= name,
            description= "empty", 
            dbxrefs=[]
            )
    print "Type 1"
    print type(rec)

    gff_File = open(gffFile, 'r')
    gffLine = []
    for line in gff_File:
        if not line.startswith('#'):
            print line
            col = line.split( )
            print col[0]
            print col[1]
            # create db if none exist
            if biodb_name != col[1] or biodb_name == "testsource":
                biodb_name = col[1]
                while True : 
                    try :
                        db = server[biodb_name] 
                        break
                    except KeyError :
                        #print ("Cannot find biodatabase with name %r making it" % biodb_name)
                        server.new_database(biodb_name)
                        server.commit()
                        #server.remove_database(biodb_name)
                        
            # controll one per gff input to global 
            if acc_id != col[0].split('.')[0]:   #test if accesseion nr set assume NC_012659.1 format	
                acc_id  = col[0].split('.')[0] 
                name    = col[0]
            
                # make seqRecord if none made so feature of rRNA can be added
                # #S   = server["refseq"].lookup(gi= gi_id).format("fasta") # add apllhabet

                S= Seq("MAAGVKQLADDRTLLDG", IUPAC.protein)

                            
                rec = SeqRecord(
                        S, 
                        id= acc_id, 
                        name= name,
                        description= description, 
                        dbxrefs=[]
                        )

                # Error duplicate key value violates unique constraint 
                # "bioentry_accession_biodatabase_id_version_key"
                # DETAIL:  Key (accession, biodatabase_id, version)=(nc_001, 57, 0) already exists.
                # meaning gi nr and accession nr must not already exist in the database

    
                # Step 2 adding gi_id
                
                #rec.annotations["gi"] = gi_id       # important for GI nummber
                    
                print "Adding SeqRecord %r to DB " % gi_id
                print biodb_name
                #db.load([rec])
                #server.commit()



            # Step 3 appending features needs 
            # TODO make a loop for 1< feature per genome 
            #rec = server[dbname].lookup(gi= gi_id)
            print "Type 2"
            print type(rec)


            f_type 	    = col[2]                # feature type NB spellig important !? !
            start 	    = int(col[3])                # start of feature (rRNA) in the sequence 
            end 	    = int(col[4])                # end of feature (rRNA) in the sequence
            #score           = col[5]                # score 
            if col[6]=='+' :
                strand  = 1
            else :
                strand = -1   # + or -, defined by 1 or -1 direction
            #frame           = col[7]                # frame
            attribute       = col[8]                #"23S ribosomal RNA"    


            product = {"product": attribute}    # soooo hard to figur out

            rec.features.append(
                    SeqFeature(
                        FeatureLocation(
                            start, 
                            end, 
                            strand= strand
                            ), 
                        type = f_type, 
                        qualifiers = product
                        )
                    )
            print "Type 3"
            print type(rec)

            #feature_location = FeatureLocation(start, stop, strand= strand)
            #rRNA_features = SeqFeature(feature_location, type = f_type, qualifiers = product) 
            #rec.features.append( rRNA_features )

    if db != "no_biodb" :
        db.load([rec])
    else :
        print acc_id
        print "not rRNA"
server.adaptor.commit()


#giID = "255534169"                        # defaul: gi id nr for NC_013062.1

acc_id = "NC_018830.1"

#gffLoader(acc_id,"RNAmmer-1.2", "r2.txt")

#gffLoader("testsource", "229599883","r2.txt") # future testing


#print server["testsource"].lookup(gi= "987654321")
#print server["RNAmmer-1.2"].lookup(gi= "255534169").features 

'''
DB = server[biodb_name]
    seq_record = DB.lookup(gi= ID) 
'''



'''

# for cleaning up gff files with seqname:
# 'gi|167234445|ref|NP_001107837. 447'
# to :
# 167234445


def newDataSorter( results, new ):
    theFile = open(results, 'r')
    newFile = open(new , 'w')
    newLine = []
    for line in theFile:
        if not line.startswith('#'):
            words = line.split( )
            tmpID= words[0].split('|')[1]
            if float(tmpID) != 0:
                words[0] = tmpID
                words.append('\n')
                #print tmpID
                new = ('\t').join(words)
                print new
                newFile.write(new)
        elif line.startswith('#'):
            newFile.write(line)

newDataSorter('results.txt', 'new.txt')
'''

'''

# for partin seq data like
# '>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655, complete genome'
def split(title):
    ID = ''
    name = ''
    desc = title
    title = desc.replace('>',"")
    words = title.split('|')
    if len(words) > 1 :
        for i in range(len(words)):
            print 'i'
            print i
            if (words[i] == 'gi') and (len(words) >= i+1) :
                ID = words[i+1]
                print ID
            if (words[i] == 'ref') and (len(words) >= i+1) :
                name = words[i+1]
                print i 
                print name 
        return ID,name,desc
    else: 
        print 'nope'


a = 'gi|167234445|ref|NP_001107837. 447'
b = '>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655, complete genome'
print split(a)
print split('fjfg|hl|k')
print split(b)
'''

''' output 
[ailaka@vestur code]$ python gffLoader.py 
Please enter user name: ******	
and password: ******
255534169	RNAmmer-1.2	rRNA	433939	436706	2858.0	+	.	23s_rRNA	

255534169	RNAmmer-1.2	rRNA	2233027	2235794	2858.0	-	.	23s_rRNA	

255534169	RNAmmer-1.2	rRNA	436856	436962	90.0	+	.	5s_rRNA	

255534169	RNAmmer-1.2	rRNA	2232771	2232877	90.0	-	.	5s_rRNA	

255534169	RNAmmer-1.2	rRNA	431637	433141	1596.8	+	.	16s_rRNA	

255534169	RNAmmer-1.2	rRNA	2236592	2238096	1596.8	-	.	16s_rRNA	

0
167234445
1
2
2
NP_001107837. 447
3
0
1
0
556503834
1
2
2
NC_000913.3
3
4
'''


