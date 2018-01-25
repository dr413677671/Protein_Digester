import re
import random
import argparse
#####################################Read Peremeters###################################################################################################################################################################################
parser = argparse.ArgumentParser(description="""2017-10-43 0:29  Ran D
# ===========================================================================
#
#               Protein_Digester
#
#  The Protein_digester script makes it easy to digest proteins.
#
#  The program is able to cope with one or more protein sequences, and digest them 
#  into small peptide sequences.The script subsequently also automatically generates 
#  the outcome data in fasta format. It can also spot unusual amino acid or unknown 
#  one, and allow user-defined cleavages. The optional enzyme in option are Trypsin,
#  Endoproteinase Lys-C, Endoproteinase Arg-C, V8 proteinase.  
#
#  As of Python >= 3.5.3, the script require some Python standard library. For 
#  users who still need to support Python 3(Anaconda),several packages are required: 
#  re(Regular expression), random and argparse, but also may support older Python 
#  versions.
#
# ===========================================================================
#
# Author:  Ran D
#
# Program Features:  
#   Optional digesting enzymes 
#   User-defined cleavages
#   Unusual or unknown amino acid identifier
#   FASTA format support
#
# ========================""",formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('files', help='input sequence files(fasta format) with '' and devided by blank',type=str, nargs='*')
parser.add_argument('-t', '--ENZYME_TYPE', help='''option of the digesting enzyme type: 1 for  Trypsin:  cuts at Lysine (Lys, K) or Arginine (Arg, R) unless the next amino acid is Proline (Pro, P), 2 for Endoproteinase Lys-C: cuts at Lysine (Lys, K) unless the next amino acid is Proline (Pro, P), 3 for  Endoproteinase Arg-C: cuts at Arginine (Arg, R) unless the next amino acid is Proline (Pro, P),  4 for V8 proteinase (Glu-C): cuts at Glutamic acid (Glu, E) unless the next amino acid is Proline (Pro, P) and default for enzyme Trypsin''', dest  =  'type', type=int)
parser.add_argument('-e', '--ERROR_TRAPPING', help='option weather to report the unusual amino acid or the unknown one(unusual amino acid as B,O,U,J,Z; unknown amino acid as X) beside peptide infor mation in the fasta format output: 1 for report 0 for not report default 0  ', dest  =  'errortrapping', type=int)
parser.add_argument('-o', '--OUTPUT', help='set the output filename', dest  =  'output', type=str)
parser.add_argument('-c', '--CLEAVAGE', help='input the number of miss cleavage of the analysis', dest  =  'cleavage', type=int)
args = parser.parse_args()
if args.files:
    print("Files List:")
    for i in args.files:
        print('%s' % (i))
else:
    print ('Load files error: Err code 1- No such files')
if args.cleavage:
    print("Cleavage:  ")
    print('%s' % (args.cleavage))
elif ((type(args.cleavage) is int)!= True):
    print ('Cleavage Error: Err code 2 -  Not a Integer')
elif (args.cleavage<0):
      print ('Cleavage Error: Err code 3 -  Smaller than zero')
errtrapping  =  0  
#####################################Read Files##################################################################################################################################################################################
def read_files(files):
    for eachfile in files:
        try:
            fasta = open('%s' % eachfile,'r')
            print ('File %s loaded successfully' % eachfile)
            fasta_sorter(fasta)
        except:
            print ('Error in loading file %s' % eachfile)
    for i in range(0,len(fa_Info),1):
        print('Sequence Information: %s \nSequence itself: %s' %(fa_Info[i],fa_Seq[i]))
#####################################Fasta Sorter##################################################################################################################################################################################
fa_num = -1
fa_Info = []#Store the seqinfo
fa_Seq = []#Store the nuc
def fasta_sorter(file):
    fa_num = -1
    info = dict()
    lines = file.readlines()
    for line in lines:
        line = line.rstrip()
        searchinfo = re.search(r'(>)([^>]+)', line, re.M| re.I)
        if searchinfo:
            fa_num += 1
            fa_Info.append(searchinfo.group(2))
            fa_Seq.append("")
        else:
            fa_Seq[fa_num] = fa_Seq[fa_num] + line
    return fa_Info, fa_Seq
#####################################Unusual_or_Unknown_identifier##################################################################################################################################################################################
def unusual_or_unknown_identifier(string):
    amino = ' '
    if (errtrapping  ==  1):  
        b = re.search(r'[BOUJZ]', string, re.I| re.M)
        if b:  
            amino = amino + '  This peptide contain unusual amino acid:  '
        b = re.search(r'B', string, re.I| re.M)
        if b:  
            amino = amino + 'B'
        b = re.search(r'O', string, re.I| re.M)
        if b:  
            amino = amino + 'O'
        b = re.search(r'U', string, re.I| re.M)
        if b:  
            amino = amino + 'U'
        b = re.search(r'Z', string,  re.I| re.M)
        if b:  
            amino = amino + 'J'
        b = re.search(r'Z', string,  re.I| re.M)
        if b:  
            amino = amino + 'Z'
        b = re.search(r'X', string, re.I| re.M)
        if b:  
            amino = amino + '  This peptide contain known amino acid:  X  '
    else:
        print  ('Default not trapping errors')
    return amino    
#####################################Digester################################################################################################################################################################################
class Digest():
    def __init__(self,fa_Info,fa_Seq):	
        self.pep_Info = []
        self.pep_Seq = []
        self.fa_Info = fa_Info
        self.fa_Seq = fa_Seq
        self.randomargument = 0
        self.output = ''
        self.count = 0
        self.amino = ''
    def cut_enzyme1(self):
        print  ('The enzyme cut KR: \n')
        for i in range(0,len(self.fa_Info),1):
            self.searchinfo = re.sub('([KR])([KR])', r'\1\2=',self.fa_Seq[i])
            self.searchinfo = re.sub('([KR])([^P])', r'\1=\2',self.searchinfo)
            self.searchinfo = self.searchinfo.split(r'=')
            if ((args.cleavage<= (len(self.searchinfo) - 1  )) & (args.cleavage >= 0)):
                self.randomargument =  random.sample(range(1, len(self.searchinfo)), args.cleavage)
                self.randomargument = sorted(self.randomargument)
                for o in self.randomargument:
                    self.searchinfo[o-1 - self.count]=self.searchinfo[o - 1 - self.count]+self.searchinfo[o - self.count]
                    del self.searchinfo[o - self.count]
                    self.count  +=1 
            else:  
                print  ('Cleavage Error: Err code 4 - the input cleavage out of length or not been input(the output is config in default cleavage 0)')
            for p in range(0,len(self.searchinfo),1):
                self.amino  =  unusual_or_unknown_identifier(self.searchinfo[p])
                self.output = self.output+ ('>%s peptide%s%s\n%s\n' % (self.fa_Info[i].rstrip('\n'), p+1, self.amino,  self.searchinfo[p]))
                self.amino  =  ''
            self.count  =  0  
        return  self.output
    def cut_enzyme2(self):
        print  ('The enzyme cut K: \n')
        for i in range(0,len(self.fa_Info),1):
            self.searchinfo = re.sub('([K])([K])', r'\1\2=',self.fa_Seq[i])
            self.searchinfo = re.sub('([K])([^P])', r'\1=\2',self.searchinfo)
            self.searchinfo = self.searchinfo.split(r'=')
            if ((args.cleavage<= (len(self.searchinfo) - 1  )) & (args.cleavage >= 0)):
                self.randomargument =  random.sample(range(1, len(self.searchinfo) ), args.cleavage)
                self.randomargument = sorted(self.randomargument)
                for o in self.randomargument:
                    self.searchinfo[o-1 - self.count]=self.searchinfo[o - 1 - self.count]+self.searchinfo[o - self.count]
                    del self.searchinfo[o - self.count]
                    self.count  +=1
            else:  
                print  ('Cleavage Error: Err code 4 - the input cleavage out of length or not been input(the output is config in default cleavage 0)')
            for p in range(0,len(self.searchinfo),1):
                self.amino  =  unusual_or_unknown_identifier(self.searchinfo[p])
                self.output = self.output+ ('>%s peptide%s%s\n%s\n' % (self.fa_Info[i].rstrip('\n'), p+1, self.amino,  self.searchinfo[p]))
                self.amino  =  ''
            self.count  =  0
        return  self.output
    def cut_enzyme3(self):
        print  ('The enzyme cut R: ')
        for i in range(0,len(self.fa_Info),1):
            self.searchinfo = re.sub('([R])([R])', r'\1\2=',self.fa_Seq[i])
            self.searchinfo = re.sub('([R])([^P])', r'\1=\2',self.searchinfo)
            self.searchinfo = self.searchinfo.split(r'=')
            if ((args.cleavage<= (len(self.searchinfo) - 1  )) & (args.cleavage >= 0)):            
                self.randomargument =  random.sample(range(1, len(self.searchinfo)), args.cleavage)
                self.randomargument = sorted(self.randomargument)
                for o in self.randomargument:
                    self.searchinfo[o-1 - self.count]=self.searchinfo[o - 1 - self.count]+self.searchinfo[o - self.count]
                    del self.searchinfo[o - self.count]
                    self.count  +=1 
            else:  
                print  ('Cleavage Error: Err code 4 - the input cleavage out of length or not been input(the output is config in default cleavage 0)')
            for p in range(0,len(self.searchinfo),1):
                self.amino  =  unusual_or_unknown_identifier(self.searchinfo[p])
                self.output = self.output+ ('>%s peptide%s%s\n%s\n' % (self.fa_Info[i].rstrip('\n'), p+1, self.amino,  self.searchinfo[p]))
                self.amino  =  ''
            self.count  = 0
        return  self.output
    def cut_enzyme4(self):
        print  ('The enzyme cut E: \n')
        for i in range(0,len(self.fa_Info),1):
            self.searchinfo = re.sub('([E])([E])', r'\1\2=',self.fa_Seq[i])
            self.searchinfo = re.sub('([E])([^P])', r'\1=\2',self.searchinfo)
            self.searchinfo = self.searchinfo.split(r'=')  
            if ((args.cleavage<= (len(self.searchinfo) - 1  )) & (args.cleavage >= 0)):            
                self.randomargument =  random.sample(range(1, len(self.searchinfo) ), args.cleavage)
                self.randomargument = sorted(self.randomargument)
                for o in self.randomargument:
                    self.searchinfo[o-1 - self.count]=self.searchinfo[o - 1 - self.count]+self.searchinfo[o - self.count]
                    del self.searchinfo[o - self.count]
                    self.count  +=1 
            else:  
                print  ('Cleavage Error: Err code 4 - the input cleavage out of length or not been input(the output is config in default cleavage 0)')
            for p in range(0,len(self.searchinfo),1):
                self.amino  =  unusual_or_unknown_identifier(self.searchinfo[p])
                self.output = self.output+ ('>%s peptide%s%s\n%s\n' % (self.fa_Info[i].rstrip('\n'), p+1, self.amino,  self.searchinfo[p]))
                self.amino  =  ''
            self.count  =  0
        return  self.output
#####################################Main##################################################################################################################################################################################
##file = open('testgenome.fasta','r')
read_files(args.files)
if (args.errortrapping  ==  1):
    errtrapping = 1
out = Digest(fa_Info, fa_Seq)
if (args.type  ==  1):  
    outputmessager = out.cut_enzyme1()  
elif (args.type  ==  2):  
    outputmessager = out.cut_enzyme2()  
elif (args.type  ==  3):  
    outputmessager = out.cut_enzyme3()  
elif (args.type  ==  4):  
    outputmessager = out.cut_enzyme4()  
else:
    print('Enzyme setup error: Err code 6 - cut in default enzyme')
if (args.output):
    file = open(args.output, 'w+')
    file.write(outputmessager)
    file.close()
else:
    print('Write file error: Err code 7 - can\'t write file')
