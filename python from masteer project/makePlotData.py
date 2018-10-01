from BioSQL import BioSeqDatabase
import sys
from operator import itemgetter

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

username = raw_input("Please enter user name: ")
password = raw_input("and password: ")

server = BioSeqDatabase.open_database(driver = "psycopg2", user = username, passwd = password, host = "dbpg-ifi-utv.uio.no", db = "rnammer")


class Genome:

    def __init__(self, ID):
        self.name = ID
        
        self.orginal    = []
        self.new        = []

    def add_org(self,data):
        self.orginal.append(data)
    
    def add_new(self,data):
        self.new.append(data)

    #def clean():
        # compare content of orginal to new
        # NB can only be don after bothe are added 

        
class Data:
    
    def __init__(self, ID, biodb_name):
        self.name = ID
        self.biodb_name = biodb_name
        
        self._5s  = []
        self._16s = []
        self._23s = []
        self.unknown = []


    def add_5s(self,start,stop):
        self._5s.append((start,stop))
    
    def add_16s(self,start,stop):
        self._16s.append((start,stop))

    def add_23s(self, start,stop):
        self._23s.append((start,stop))
    
    def add_unknown(self, start,stop):
        self.unknown.append((start,stop))

    def __repr__(self):
        return repr(self.biodb_name)

#
# make for all
#

class Stats:
    
    
    def __init__(self, biodb_name_one,biodb_name_two):
        self.biodb_name_one = biodb_name_one 
        self.biodb_name_two = biodb_name_two
        #self.ID         = ID
        self.genomes    = []
        #self.data       = []

    def add_genome(self,ID):
        g= Genome(ID)
        g.add_org(self.add_data(self.biodb_name_one, ID ))
        g.add_new(self.add_data(self.biodb_name_two, ID ))
        
        self.genomes.append(g)


    def add_data(self,db,ID):
        db = server[db]                		
        seq_record = db.lookup(accession= ID)

        features = seq_record.features                      # tmp to see what is happening
        rRNAs = [x for x in features if x.type == "rRNA"]   # list of <class 'Bio.SeqFeature.SeqFeature

        List = []
        data = Data(ID,db)
        
        import re
        for part in rRNAs:

            tmp2 = str(part.location)
            if "join" in tmp2:
                print "Join rRna skiped for ID : {}".format(ID)
                print tmp2
            else :
                tmp2 = re.split(r'[()\[\]:]', tmp2)
                tmp2 = filter(None,tmp2)

                att = (str(part.qualifiers["product"]))#.split()[0][2:])
                
                tmp2.append(att)#.split('_')[0])
        
                List.append(tmp2)

        #return List

    #def correct(List):
        for l in List:
            #print l[0]
            l[0] = float(l[0])
            l[1] = float(l[1])
            
            if l[2] == '-':
                tmp     = l[0]
                l[0]    = l[1]
                l[1]    = tmp
                l[2]    = '+'
        #return List

    #def attr(List):
        '''
        L1 = []
        L2 = []
        L3 = []
        U = []
        '''

        for l in List:
            #print l[3]
            if '5' in l[3]:
                #print "match 5s"
                data.add_5s(l[0],l[1])
                #L1.append(l)
            elif '16' in l[3]:
                data.add_16s(l[0],l[1])
                #L2.append(l)
            elif '23' in l[3]:
                data.add_23s(l[0],l[1])
                #L3.append(l)
            else :
                print "No attribute match ({}) to 5s,16s or 23s".format(l[3])
                data.add_unknown(l[0],l[1])
                #U.append(l)
        
        #print "5s before sorter------------------------------- "
        #print data._5s
        #print data._5s[0][0]

        data._5s.sort( key=lambda x: x[0] )
        #print "5s after sorter-------------------------------- "
        #print data._5s

        data._16s.sort( key=lambda x: x[0] )
        data._23s.sort( key=lambda x: x[0] )
        
        return data

    def makePlotData(self):
        print "PLOTTING"
        plot_data = '_'.join([self.biodb_name_one, "plot_data.txt"])

        f = open(plot_data,'w')

        diff_5s_start       = [] 
        diff_16s_start      = []
        diff_23s_start      = []

        diff_5s_stop        = [] 
        diff_16s_stop       = []
        diff_23s_stop       = []


        for g in self.genomes:
            print g.name
            old = g.orginal[0]
            new = g.new[0]


            if len(old._5s) == len(new._5s):
                for i in calcDiff(old._5s,new._5s):
                    diff_5s_start.append(i[0])
                    diff_5s_stop.append(i[1])
            #elif len(old_5s) < len(new_5s):
                # search trough olds Unkonwn to get missing
            if len(old._16s) == len(new._16s):
                for i in calcDiff(old._16s,new._16s):
                    diff_16s_start.append(i[0])
                    diff_16s_stop.append(i[1])
        
            if len(old._23s) == len(new._23s):
                for i in calcDiff(old._23s,new._23s):
                    diff_23s_start.append(i[0])
                    diff_23s_stop.append(i[1])
            print len(diff_5s_start)
        
        f.write("# 5s start")
        f.write("\n")

        for it in diff_5s_start:
            f.write("%s " %it)
        
        f.write("\n")
        f.write("# 5S stop")
        f.write("\n")

        for it in diff_5s_stop:
            f.write("%s " %it)
        
        f.write("\n")
        f.write("# 16s start")
        f.write("\n")

        for it in diff_16s_start:
            f.write("%s " %it)
        
        f.write("\n")
        f.write("# 16S stop")
        f.write("\n")

        for it in diff_16s_stop:
            f.write("%s " %it)

        f.write("\n")
        f.write("# 23s start")
        f.write("\n")

        for it in diff_23s_start:
            f.write("%s " %it)

        f.write("\n")
        f.write("# 23S stop")
        f.write("\n")

        for it in diff_23s_stop:
            f.write("%s " %it)


        print "-------------Diff------------------"
        print "16s"
        print diff_16s_start
        print diff_16s_stop




# sjekk sorter and writ into data poibt 
# method for diffrence where to put
# make reader of all 
# cleaner in stat
#

def calcDiff(ListA,ListB):
    return [(x[0]-y[0], x[1]-y[1]) for x,y in zip(ListA, ListB)]


#ex =["NC_011750","NC_004603","NC_003197"]# , "NC_004668","NC_000913"]
  
#t =  ["NC_014644","NC_003197"]  # problem child
#biodb1 = "refseq"
#biodb2 = "RNAmmer-1.2"
#ID= "NC_002967"


def stk(DBs):
    for biodb_org in DBs:
        biodb_new = '_'.join(["RNAmmer-1.2",biodb_org])
        
        S   = Stats(biodb_org,biodb_new)
        n   = 0
        db  = server[biodb_new]
        print "This database contains %i records" % len(db)
        count = 1
        #for i in t: 
        #Taking acc. id from new because if its not in that one its nothing to compare 
        for key, record in db.iteritems():
            print "Orginal id %s" % record.id
            S.add_genome(record.id)
            #print i
            #S.add_genome(i)
            print S.genomes[n].name
            print len(S.genomes[n].orginal[0]._16s)
            print "of old and new seq pred  "
            print len(S.genomes[n].new[0]._16s)
            n = n+1
        S.makePlotData()

biodb_names     = ["Actinobacteria","Bacteroidetes","Chlamydiae","Cyanobacteria","Deinococcus-Thermus","Firmicutes","Spirochaetes","Tenericutes"]

singles         = ["Aquificae","Chlorobi","Chloroflexi","Dictyoglomi","Fusobacteria","Nitrospirae","Planctomycetes","Synergistetes","Thermotogae"]

Proteobacteria  = ["Alphaproteobacteria","Betaproteobacteria","Gammaproteobacteria","Deltaproteobacteria","Epsilonproteobacteria"]

import time
start_time = time.time()
stk(Proteobacteria)
print("--- %s seconds ---" % (time.time() - start_time))


