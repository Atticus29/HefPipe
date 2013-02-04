##Created by Mark Fisher (email: mark.aaron.fisher@gmail.com)
import re
import os
import csv
import pickle
import copy
import math
import itertools
from pyper import *
#os.chdir("/Users/markfisher/Desktop")

#opens spreadsheet and assigns it to rows
def readCsv():
        filename=raw_input("Name of csv file you want to read")
        csv_reader=csv.reader(open(filename, "rU"))
        try:
                rows=[row for row in csv_reader]
        except:
                rows=[]
        return rows

def readCsv_input(filesname):
        csv_reader=csv.reader(open(filesname, "rU"))
        try:
                rows=[row for row in csv_reader]
        except:
                rows=[]
        return rows

def readCsv_input_safer(filename):
        csv_reader=csv.reader(open(filename, 'rU'))
        #print "csv_reader", csv_reader
        data=[]
        for row in csv_reader:
                if row:
                        data.append(row)
        return data
                        

## {{{ http://code.activestate.com/recipes/410687/ (r4)
def transposed(lists):
   if not lists: return []
   return map(lambda *row: list(row), *lists)

def transposed2(lists, defval=0):
   if not lists: return []
   return map(lambda *row: [elem or defval for elem in row], *lists)
## end of http://code.activestate.com/recipes/410687/ }}}


def transposeCsv(rows): #assumes that all of the rows are the same length
        temprow=[]
        cols=[[]]
        for i in range (0, len(rows[0])):
                for j in range (0, len(rows)):
                        temprow.append(rows[j][i])
                cols.append(temprow)
                temprow=[]
        #print cols
        return cols

#saves a table of tables as a csv file
def saveCsv(cols, filename, destination_directory):
        savefile=destination_directory
        savefile+=filename
        ofile=open(savefile, "wb")
        writer=csv.writer(ofile, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONE)
        for row in cols:
                #print row
                #print cols.index(row)
                writer.writerow(row) #if it throws an escapechar error, inspect the spreadsheet for weird characters such as ,;*-
                #ofile.writelines(row)
        ofile.close()

def saveCsv_v2(cols, filename, destination_directory):
        savefile=destination_directory
        savefile+=filename
        ofile=open(savefile, "wb")
        writer=csv.writer(ofile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for row in cols:
                #rint row
                #print cols.index(row)
                writer.writerow(row) #if it throws an escapechar error, inspect the spreadsheet for weird characters such as ,;*-
                #ofile.writelines(row)
        ofile.close()


def Diploid_genotype_table_to_allele_table(table): #assumes a header row with ID in first column and locus name in each following column (only once); missing data as ***s
        #make the first column ID
        newTable=[]
        newTable.append([table[0][0]])

        #now add to that header row, making sure to repeat
        for x in range(1,len(table[0])):
                       newTable[0].append(table[0][x])
                       newTable[0].append(table[0][x])

        #now go through each row of table beyond the first, adding a new row and splitting the genotypes
        for x in range (1, len(table)):
                newTable.append([table[x][0]])
                for y in range (1, len(table[x])):
                        newTable[x]+=split(table[x][y])
        #print newTable
        return newTable

def Table_to_DictOfLists(table): #assumes each header row with ID in first col and each other cell in top row populated with locus name (ie. doubled), missing data contains **, allele values are at least two digits
        samples={}
        table_trans=transposed(table)
        #print table_trans[0]
        #print len(table_trans[0])
        #print table
        locusNameList=[[table[0][x]] for x in range (1,len(table[0]),2)]
        #print locusNameList
        for x in range(1, int(table_trans[0][len(table_trans[0])-1])+1):
                if (str(x) in table_trans[0]):
                        #print x
                        samples[x]=copy.deepcopy(locusNameList)
                #print samples[x]
        #print samples
        #print samples[1]
        for y in samples:#for each individual
                #print 'y',y
                #print samples[y]
                #print table_trans[0].index(str(y))
                #go through the loci:
                for z in range(1, len(samples[1])+1):
                        #print 'y',y
                        #print 'z',z
                        #print type(samples[y][z-1])
                        #print samples[y][z-1]
                        samples[y][z-1].append(table[table_trans[0].index(str(y))][table[0].index(samples[y][z-1][0])])
                        samples[y][z-1].append(table[table_trans[0].index(str(y))][table[0].index(samples[y][z-1][0])+1])
                        #print samples [y][z-1]
        return samples                                       
        

def DictOfLists_toTable(dictOfLists, pipeline_directory, rejected_samples_address):
        #print "does it happen before it gets passed to DictOfLists_toTable?", dictOfLists
        #print "dictOfLists['371']", dictOfLists['371']

        #check for prolematic locus still hanging around
        #for (key, lst) in dictOfLists.items():
         #       for locus in lst:
          #              if len(locus)!=3:
                                #print "the key %s locus %s" %(key, locus)
                
        address_of_rejected_samples=rejected_samples_address
        exceptions=readCsv_input(address_of_rejected_samples)
        data=[['ID']]
        #if ([str(int(x))] not in exceptions):
        #populate the top row with locus names
        #for x in dictOfLists:
        for i in range (0, len(dictOfLists[dictOfLists.keys()[0]])):
                data[0].append(dictOfLists[dictOfLists.keys()[0]][i][0])
                data[0].append(dictOfLists[dictOfLists.keys()[0]][i][0]) #twice because the locus name should appear once for each allele

        #go through all of the other rows
        for x in range (0, len(dictOfLists)):
                #print "x", x
                data.append([])
                data[len(data)-1].append(dictOfLists.keys()[x])# make sure these samples match their genotypes
                for y in range(0,len(dictOfLists[dictOfLists.keys()[0]])):
                        #print "y", y
                        if (dictOfLists.keys()[x] not in exceptions and len(dictOfLists[dictOfLists.keys()[x]][y])==3):
                                #print "this happens at x", x, " and y", y
                                data[len(data)-1].append(dictOfLists[dictOfLists.keys()[x]][y][1])#just changed this from x+1 to x for indexing errors after generalizing. 6.24.2012
                                data[len(data)-1].append(dictOfLists[dictOfLists.keys()[x]][y][2])#see lijne above
        data=remove_spaces_from_list_recursive(data)
        saveCsv(data, 'final_output.csv', pipeline_directory)
        return data

def GenePop_ForLDNE (dictOfLists, pipeline_directory, rejected_samples_address): #this function works like GenePop_Output, except all of the alleles will be the same length (ie. 84 83 genotype will be 084083 if there's a 100 allele at some other locus!)
        #assume no exceptions at this point (assumes a fully fleshed out dataset), also assumes that mis-scores or failures are '**'
        data=GenePop_Output(dictOfLists, pipeline_directory, rejected_samples_address)
        #print data
        exceptions=[]
        #dictOfLists=clean_trouble_loci(dictOfLists, exceptions)
        global_diploid_genotype_maxlength=Max_length_of_any_genotype(dictOfLists)

        #now go through every value in data and if it's not ml insert a 0 at the beginning of it
        for x in range(0, len(data)):
                for y in range (2,len(data[0])):
                        #print data[x][y]
                        #print len(data[x][y])
                        if (len(data[x][y])<global_diploid_genotype_maxlength):
                                data[x][y]=insert(data[x][y],'0',0)
                                #print "did it"
        #now compress the data
        for x in range (0, len(data)):
                for y in range (2, len(data[0])-1,2):
                        data[x][y]=insert(data[x][y],data[x][y+1],3)
        #print data
        newdata=[]
        for x in range(0, len(data)):
                newdata.append(data[x][0:1]+['_'])
                for y in range(2, len(data[0])-1,2):
                        newdata[x].append(data[x][y])

        for x in range (0, len(newdata)):
                for y in range (0, len(newdata[0])):
                        if ('*' in newdata[x][y]):
                                newdata[x][y]='000000'
        saveCsv(newdata,str(dictOfLists[dictOfLists.keys()[1]][1])+'_LDNE_ready.txt', pipeline_directory)
        print [dictOfLists[dictOfLists.keys()[1]][x][0] for x in range(0,len(dictOfLists[dictOfLists.keys()[1]]))]
        
def GenePop_Output(dictOfLists, pipeline_directory, rejected_samples_address):
        exceptions=readCsv_input(rejected_samples_address) # a list of lists of strings
        #print "exceptions", exceptions, type(exceptions)
        #print "dictOfLists as it gets passed to GenePop_Output", dictOfLists
        keys_list=dictOfLists.keys()
        keys_as_numbers=[]
        for x in keys_list:
                if x != 'Blank' and x!='blank':
                        keys_as_numbers.append(int(x))
        #print "key_list", keys_list
        keys_as_numbers_sorted=sorted(keys_as_numbers)
        #print "now sorted", keys_as_numbers_sorted

                ##left off here!
                
        #for when you want exceptions to be input by the user rather than taken from a file, the syntax is in the next few lines below
        data=[]
        #print dictOfLists[3]
        dictOfLists=clean_trouble_loci(dictOfLists, exceptions)
        #print "after cleaning trouble loci", dictOfLists
        final_output_address=pipeline_directory
        saveFile= open(final_output_address+'final_output_GenePop.txt', 'w')
        saveFile.write(raw_input("Population name?"))
        saveFile.write('\n')
        #write the locus names
        for n in range(0, len(dictOfLists[dictOfLists.keys()[1]])):#obviously, this should be 1 or could be any entry in the dictionary, but I know that there are different numbers of loci being scored right now, so I picked a later one to be safe
                if (len(dictOfLists[dictOfLists.keys()[1]][n])>1):                  
                    saveFile.write(dictOfLists[dictOfLists.keys()[1]][n][0].replace(' ', '_'))
                    saveFile.write('\n')
        saveFile.write('POP'+'\n')

        #write the rest of the file if the dictionary keys are sorted numbers
#un-indented this for loop
        for a in keys_as_numbers_sorted:
                #print "a", a
                #print "type(a)", type(a)
                data.append([])
                data[len(data)-1].append(str(a)) #inserts the sample name
                if ([str(a)] not in exceptions and str(a) not in exceptions):
                        saveFile.write(str(a))
                        saveFile.write(', ')#6.25.2012 added a space here
                data[len(data)-1].append(', ') # inserts a comma after the sample name (BUT HOW TO DEAL WITH THE SPACE??)
                #saveFile.write(',')
                for y in range(0, len(dictOfLists[dictOfLists.keys()[1]])):# same as len(keeplist)
                        #print "y", y
                        if (str(a) not in exceptions and [str(a)] not in exceptions and len(dictOfLists[str(a)])>1):
                                #print "yth locus", dictOfLists[a][y]
                                if (len(dictOfLists[str(a)][y])>1):
                                        data[len(data)-1].append(dictOfLists[str(a)][y][1]) #how to get rid of the inevitable comma between this and what follows it
                                        #print "dictOfLists [a][y]", dictOfLists[a][y]
                                        #print 'dictOfLists[str(a)][y][1] and [2]', dictOfLists[str(a)][y][1], dictOfLists[str(a)][y][2]
                                        if(dictOfLists[str(a)][y][1]=='**' or dictOfLists[str(a)][y][1]=='0**' and dictOfLists[str(a)][y][2]!='**'):
                                                #print "this should be a homozygote:", str(dictOfLists[str(a)][y][1]+dictOfLists[str(a)][y][2])
                                                replace_zeros1=(str(dictOfLists[str(a)][y][2]+dictOfLists[str(a)][y][2])+' ').replace('*','0')
                                                saveFile.write(replace_zeros1)

                                        else:
                                                if(dictOfLists[str(a)][y][2]=='**' or dictOfLists[str(a)][y][2]=='0**' and dictOfLists[str(a)][y][1]!='**'):
                                                        #print "this should be a homozygote:", str(dictOfLists[str(a)][y][1]+dictOfLists[str(a)][y][2])
                                                        #print "this locus is in", dictOfLists[str(a)][y]
                                                        replace_zeros=(str(dictOfLists[str(a)][y][1]+dictOfLists[str(a)][y][1])+' ').replace('*','0')
                                                        saveFile.write(replace_zeros)

                                                else: #this still might be problematic
                                                        if (str(dictOfLists[str(a)][y][1]+dictOfLists[str(a)][y][2])=='****' or (dictOfLists[str(a)][y][1]=='**' or dictOfLists[str(a)][y][2]=='**' or dictOfLists[str(a)][y][1]=='0**' or dictOfLists[str(a)][y][2]=='0**') and len(str(dictOfLists[str(a)][y][1]+dictOfLists[str(a)][y][2]))!=5):
                                                                for q in range(0, Max_length_of_diploid_genotype(dictOfLists, y, exceptions)):
                                                                        saveFile.write('0')
                                                                saveFile.write(' ')
                                                                #print "this locus is in", dictOfLists[x][y]
                                                        else:
                                                                if (len(str(dictOfLists[str(a)][y][1]+dictOfLists[str(a)][y][2]))==5):
                                                                        saveFile.write('000000 ')
                                                                        #print "this locus is in", dictOfLists[x][y]
                                                                        #print "problem with sample ", x+1, " at locus", dictOfLists[x+1][y][0], " (the ", y, "th locus)"
                                                                else:
                                                                        saveFile.write(str(dictOfLists[str(a)][y][1]+dictOfLists[str(a)][y][2])+' ')   #will have to keep track of the maximum length of this, run another for loop and change all values missing this max length with 0s equal in number to max length (warning-loci with one 2-digit allele and another 3-digit allele will basically be lost-I think that they will be reported as uninformative loci, though)                     
                                        data[len(data)-1].append(dictOfLists[str(a)][y][2])
                                #saveFile.write(str(dictOfLists[a+1][y][1]))
                saveFile.write('\n')
        saveFile.close()

        #get rid of blank lines
        rid_of_blanks=open(pipeline_directory+"final_output_GenePop.txt").readlines()
        rid_of_blanks=filter(lambda x: not x.isspace(), rid_of_blanks)
        save_again=open(pipeline_directory+"final_output_GenePop.txt", "w")
        for line in rid_of_blanks:
                save_again.write(line)
        return keys_as_numbers_sorted

def clean_trouble_loci(dictOfLists, exceptions): #takes loci like Sway and i-136, which have some alleles that are 2 digits long (82) and some that are 3 digits long (100)
        #print "dictOfLists", dictOfLists
        for x in dictOfLists:
                #print "x", x
                #print dictOfLists[x]
                #print "len of dictOfLists[x]", len(dictOfLists[x])
                if (x not in exceptions and len(dictOfLists[x])>1):# (modified for Table_to_DictOfLists function 11.28.2011)(modified again, along with much of the script, for generalization in 6.2012)
                        #print "not in exceptions happens"
                        ml=0
                        for y in range (0, len(dictOfLists[dictOfLists.keys()[0]])):
                                #print "y", y
                                #print dictOfLists[x][y][1]
                                #print dictOfLists[x][y][2]
                                if (len(dictOfLists[x][y])>1): #I think I had to add this to bypass loci that are in keeplist that don't exist in this particular analysis??
                                        #print "x", x, "y", y, "dictOfLists[x][y]", dictOfLists[x][y]
                                        if(len(dictOfLists[x][y][1])+len(dictOfLists[x][y][2])<Max_length_of_diploid_genotype(dictOfLists,y,exceptions)):
                                                if (len(dictOfLists[x][y][1])<len(dictOfLists[x][y][2])):#huge big o here
                                                        #print "1<2 insert happens"
                                                        dictOfLists[x][y][1]=insert(dictOfLists[x][y][1],'0',0)
                                                else:
                                                        if (len(dictOfLists[x][y][2])<len(dictOfLists[x][y][1])):
                                                                #print "2<1 insert happens"
                                                                dictOfLists[x][y][2]=insert(dictOfLists[x][y][2],'0',0)
                                                        else:
                                                                dictOfLists[x][y][2]=insert(dictOfLists[x][y][2],'0',0)
                                                                dictOfLists[x][y][1]=insert(dictOfLists[x][y][1],'0',0)
                                                        
        return dictOfLists
                                                

def Max_length_of_diploid_genotype(dictOfLists, y, exceptions):
        ml=0
        #print "len of DictOfLists", len(dictOfLists)
        #print "the y that got passed to this Max_length_of_diploid_genotype function", y
        for x in dictOfLists:
                #print "x+1", x+1
                #print "x", x
                if x not in exceptions: #this might now be wrong (modified for Table_to_DictOfLists function 11.28.2011)
                        #print "x+1", x+1, "y", y
                        #print dictOfLists[x+1][y]
                        if(len(dictOfLists[x][y])>1):
                                #print "dictOfLists[x][y] inside Max_length_of_diploid_genotype function", dictOfLists[x][y]
                                #print "dictOfLists[x][y][1]", dictOfLists[x][y][1]
                                #print "dictOfLists[x][y][2]", dictOfLists[x][y][2]
                                if('*' not in dictOfLists[x][y][1] and '*' not in dictOfLists[x][y][2]):
                                        if (len(dictOfLists[x][y][1]+dictOfLists[x][y][2])>ml):
                                                ml=len(dictOfLists[x][y][1]+dictOfLists[x][y][2])
                                                #print "sample", x, "locus", dictOfLists[x][y], "max length", ml
        return ml

def Max_length_of_any_genotype(dictOfLists):#assumes no exceptions-same assumptions as GenePop for LDNE function
        ml=0
        for x in dictOfLists:
                for y in range (0, len(dictOfLists[dictOfLists.keys()[1]])):
                        if (len(dictOfLists[x][y][1]) >ml):
                                ml=len(dictOfLists[x][y][1])
                        if (len(dictOfLists[x][y][2])>ml):
                            ml=len(dictOfLists[x][y][2])
        #print "max length of any genotype is ",ml
        return ml
        
                
def insert(original, new, pos):
#Inserts new inside original at pos for strings
        return original[:pos] + new + original[pos:]

def assign_zygosity(string):#for the convert to g2 script
        if (len(string)%2==0 and not special_match(string)): #this would mean the number of digits is even in the entry, suggesting that it's not a weird entry
                if(string[:len(string)/2]==string[len(string)/2:]):
                        return int(0)
                else:
                        return int(1)
        else:
                #print "this happens"
                if (special_match(string)):
                        return int(-99)
                else:
                        if(string[(len(string)-1)/2]=='0'):
                                tempstring=string[:((len(string)-1)/2)]+string[((len(string)-1)/2+1):]#gets rid of middle digit?? 4_8_2012
                                #print"tempstring", tempstring
                                return assign_zygosity(tempstring)
                        else:
                                return int(1)

def split(locus): #gets passed a string of x digits representing a diploid genotype, some with five digits
        choiceA=[locus[0:len(locus)/2+1],locus[len(locus)/2+1:len(locus)]]
        choiceB=[locus[0:len(locus)/2],locus[len(locus)/2:len(locus)]]
        if('*' in choiceA[0] or '*' in choiceA[1]):
                return choiceB
        if(locus=='0'):
                return [0,0]
        if (math.fabs(int(choiceA[0])-int(choiceA[1]))<math.fabs(int(choiceB[0])-int(choiceB[1]))):
                return choiceA
        else:
                return choiceB

def splitUnderscore(locus): #gets passed a string of x digits representing a diploid genotype, some with five digits
        choiceA=locus[0:len(locus)/2+1],"_",locus[len(locus)/2+1:len(locus)]
        choiceB=locus[0:len(locus)/2],"_",locus[len(locus)/2:len(locus)]
        if('*' in choiceA[0] or '*' in choiceA[1]):
                return choiceB
        if(locus=='0'):
                return "0_0"
        if (math.fabs(int(choiceA[0])-int(choiceA[1]))<math.fabs(int(choiceB[0])-int(choiceB[1]))):
                return choiceA
        else:
                return choiceB

def rhh_output(data, pipeline_directory):
        #print "data being passed to rhh", data
        saveFile=open(pipeline_directory+'rhh_ready.txt', 'w')
        for x in range (0, len(data)):
                if (len(data[x])>0):
                        #print "data [",x,"]", data[x]
                        saveFile.write(data[x][0]+" ")
                        for z in range(0, len(data[x][1])):
                                        #print "data[x][y][z]", data[x][y][z]
                                        #print "type", type(data[x][y][z])
                                        saveFile.write(data[x][1][z][0]+ " ")
                                        if (len(data[x][1][z])>1):
                                                saveFile.write(data[x][1][z][1]+" ")
                        saveFile.write('\n')
        saveFile.close()
        #print "You may want to make sure that all of your samples have the same number of genotypes in rhh_ready.txt in the pipeline folder (the rhh script will also throw you an error, but just in case you want things to run smoothly)"
        return data

def list_of_lists_to_txt (data, desktop):
        name=raw_input("Name of the file you want to save (default used to be data_as_text_file)...will automatically save to desktop as a .txt file, so those parts are unnecessary. To remind you, this will be saving a list of lists (i.e. a data table) into a text file (perhaps to be passed on to a new function in a second.")
        print "\nSo the file will be:"+desktop+name+".txt" 
        saveFile=open(desktop+name+'.txt', 'w')
        for x in range (0, len(data)):
                for y in range(0, len(data[x])):
                        #print "data [",x,"]", data[x]
                        saveFile.write(data[x][y]+" ")
                        
                saveFile.write('\n')
        saveFile.close()
        print "You may want to make sure that all of your samples have the same number of genotypes (the rhh script will also throw you an error, but just in case you want things to run smoothly)"
        return data

def txt_to_csv(filename):
        #print "filename", filename
        #print type(filename)
        inputFile=open(filename).readlines()
        #print "inputFile", inputFile
        inputFileList= [ line.strip().split('\t') for line in open(filename)]
        #print "inputFileList", inputFileList
        writer=csv.writer(file(filename+'_csv_converted.csv', 'wb'))
        writer.writerows(inputFileList)
        return (writer)
def txt_to_csv_space_split_version(filename):
        #print "filename for txt_to_csv_space_split_version", filename
        
        inputFileList= [ line.strip().split() for line in open(filename, 'r')]
        #print "inputFileList", inputFileList
        writer=csv.writer(file(filename+'_csv_converted.csv', 'wb'))
        writer.writerows(inputFileList)
        return (writer)
def remove_two_word_locus_names(filename):#should also get rid of any non-number entries beyond the first column, actually (e.g. "switches")
        data=readCsv_input(filename)
        #print "data from remove_two_word_module", data
        remove=[]
        for x in range (0, len(data)):
                #print "len(data[x])", len(data[x])
                for y in range (1, len(data[x])):
                        if len(data[x])>1:
                                #print "x", x, "\ny", y, "\n data[x][y]", data[x][y]
                                if re.match("^[A-Za-z]*$", data[x][y]):
                                        #print "this happens"
                                        #print "address of", data[x][y], " gets stored"
                                        remove.append([x,y])
        #print "remove", remove
        for z in range (0, len(remove)):
                #print "z", z
                #print "len of remove", len(remove)
                del data[remove[len(remove)-1-z][0]][remove[len(remove)-1-z][1]]
        #print "after editing", data
        name=list_of_lists_to_txt_sensitivity_specific (data, filename)
        print "name: ", name
        #name=raw_input("Just paste the filename that just got printed here")
        txt_to_csv_space_split_version(name)

def list_of_lists_to_txt_sensitivity_specific (data, filename):
        name=filename+'data_as_text_file.txt'
        print "\nSo the file will be:"+name #take this out
        saveFile=open(name, 'w')
        for x in range (0, len(data)):
                for y in range(0, len(data[x])):
                        #print "data [",x,"]", data[x]
                        saveFile.write(data[x][y]+" ")
                        saveFile.write('\n')
        saveFile.close()
        #print "You may want to make sure that all of your samples have the same number of genotypes (the rhh script will also throw you an error, but just in case you want things to run smoothly)"
        return name

def testDict(dictOfLists, rejected_samples_address):
        address_of_rejected_samples=rejected_samples_address
        exceptions=readCsv_input(address_of_rejected_samples)
        data=[['ID']]
        #if ([str(int(x))] not in exceptions):
        #populate the top row with locus names
        for x in dictOfLists:
                print x
        return data
                
def special_match(strg, search=re.compile(r'[^0.]').search): #http://stackoverflow.com/questions/1323364/in-python-how-to-check-if-a-string-only-contains-certain-characters
        return not bool(search(strg))

def allele_report_pipeline (allele_reports_address, keeplist_address, monolist_address, desktop, pipeline_directory, rejected_samples_address):
    exceptions=readCsv_input(rejected_samples_address)
    #note: you will have to run this more than once. The second time will be after you have dealt with all of the issues related to missing.txt
    alleleList=readCsv_input(allele_reports_address)
    present_cumulative=[]
    keeplist=readCsv_input(keeplist_address)
    #print "keeplist",keeplist
    monolist=readCsv_input(monolist_address)
    #print "monolist", monolist
    renamed_files=[]


    #go through all of the allele reports on the allelereports list and do the following things:    
    for a in range(0, len(alleleList)): #len(alleleList)
        #print "a", a
        #open the file
        file=desktop
        file+=(alleleList[a][0])
        print "file", file
        data=readCsv_input(file)
        
        #get rid of the top rows with the unimportant information
        #print "data before cleanup", data
        cleaned=cleanup(data)
        #print "data after cleanup", cleaned

        #fill the top row with names corresponding to loci
        filled=fill_header(cleaned)
        #print "data after fill_header", filled

        #print "The file",alleleList[a][0], "needs to have the following loci fixed:"
        #remove loci that have three or more columns ----------------This is now also the spot/function to check for items on keeplist that are being removed
        removed_bad=remove_bad_loci(filled, keeplist)
        #print "after bad removed", removed_bad

        #keeps track of every locus that has been processed
        for x in range (0, len(removed_bad[0])): #this assumes that the locus names are in the first row, which they should be at this point
            present_cumulative.append(removed_bad[0][x])
        #print "present_cumulative", present_cumulative
        #gets rid of loci still in the samples that we know apriori (from monolist) are monomorphic
        rid=remove_nonscores(removed_bad, monolist)

        #renames the ID with numbers instead of their long crazy names
        rid_trans=transposed(rid)
        #print "rid_trans", rid_trans
        delimiter=raw_input("If your sample names have a number in the front of them separated by a delimiter, what is that delimiter? Otherwise, press enter.")
        if len(delimiter) !=0:
            #print "before renaming row", rid_trans[1]
            rid_trans[1]=rename_row(rid_trans[1],delimiter)
            #print "after renaming row", rid_trans[1]
        rid_back=transposed(rid_trans)
        
        filename='edited_'
        filename+=alleleList[a][0]
        renamed_files.append([filename])
        saveCsv(rid_back,filename, pipeline_directory)

    #cross references the loci that have been scored with a list of loci that SHOULD have been scored and reports any on the keeplist that were missed (usually means some were scored as triploids, etc.)    
    most_missingest(present_cumulative, keeplist, pipeline_directory)
    #print renamed_files
    saveCsv(renamed_files, 'edited_file_list', pipeline_directory)

    #check that all of the samples that are present in one file are present in the others
    unique_samples_in_edited_files(pipeline_directory+'edited_file_list', pipeline_directory, exceptions)

def allele_report_pipeline_directory_version (allele_reports_address, keeplist_address, monolist_address, pipeline_directory, rejected_samples_address):
    exceptions=readCsv_input(rejected_samples_address)
    #note: you will have to run this more than once. The second time will be after you have dealt with all of the issues related to missing.txt
    #import os
    #import csv
    ##mport allele_report_manip
    #import transpose
    alleleList=[x for x in os.listdir(allele_reports_address) if '.csv' in x]
    #print "alleleList", alleleList
    present_cumulative=[]
    keeplist=readCsv_input(keeplist_address)
    #print "keeplist",keeplist
    monolist=readCsv_input(monolist_address)
    #print "monolist", monolist
    renamed_files=[]


    #go through all of the allele reports on the allelereports list and do the following things:    
    for a in alleleList: #len(alleleList)
        #print "a", a
        #open the file
        file=allele_reports_address
        file+=a
        print "file", file
        data=readCsv_input(file)
        
        #get rid of the top rows with the unimportant information
        #print "data before cleanup", data
        cleaned=cleanup(data)
        #print "data after cleanup", cleaned

        #fill the top row with names corresponding to loci
        filled=fill_header(cleaned)
        #print "data after fill_header", filled

        #print "The file",alleleList[a][0], "needs to have the following loci fixed:"
        #remove loci that have three or more columns ----------------This is now also the spot/function to check for items on keeplist that are being removed
        removed_bad=remove_bad_loci(filled, keeplist)
        #print "after bad removed", removed_bad

        #keeps track of every locus that has been processed
        for x in range (0, len(removed_bad[0])): #this assumes that the locus names are in the first row, which they should be at this point
            present_cumulative.append(removed_bad[0][x])
        #print "present_cumulative", present_cumulative
        #gets rid of loci still in the samples that we know apriori (from monolist) are monomorphic
        rid=remove_nonscores(removed_bad, monolist)
        #print "rid", rid

        #renames the ID with numbers instead of their long crazy names
        rid_trans=transposed(rid)
        #print "rid_trans", rid_trans
        delimiter=raw_input("If your sample names have a number in the front of them separated by a delimiter, what is that delimiter? Otherwise, press enter.")
        if len(delimiter) !=0:
            #print "before renaming row", rid_trans[1]
            rid_trans[1]=rename_row(rid_trans[1],delimiter)
            #print "after renaming row", rid_trans[1]
        rid_back=transposed(rid_trans)
        #print "rid_back", rid_back
        filename='edited_'
        filename+=a
        renamed_files.append([filename])
        saveCsv(rid_back,filename, pipeline_directory)

    #cross references the loci that have been scored with a list of loci that SHOULD have been scored and reports any on the keeplist that were missed (usually means some were scored as triploids, etc.)    
    most_missingest(present_cumulative, keeplist, pipeline_directory)
    #print renamed_files
    saveCsv(renamed_files, 'edited_file_list', pipeline_directory)

    #check that all of the samples that are present in one file are present in the others
    unique_samples_in_edited_files(pipeline_directory+'edited_file_list', pipeline_directory, exceptions)


def unique_samples_in_edited_files(edited_files_address, pipeline_directory, exceptions):
    edited_file_list=readCsv_input(edited_files_address)
    #print "edited_file_list", edited_file_list
    entire_list=[]

    #make a list of every sample name that occurs
    for x in range (0, len(edited_file_list)):
            current_file=readCsv_input(pipeline_directory+edited_file_list[x][0])
            #print "file name", edited_file_list[x][0]
            #print "current_file", current_file
            for y in range (1, len(current_file)):
                    if current_file[y][1] not in entire_list:
                            entire_list.append(current_file[y][1])
    #print "entire_list", entire_list

    for x in range (0, len(edited_file_list)):
            current_file=readCsv_input(pipeline_directory+edited_file_list[x][0])
            print "file name", edited_file_list[x][0]
            list_of_samples_in_edited_file=[current_file[y][1] for y in range(1, len(current_file))]
            #print "list_of samples_in_edited_file", list_of_samples_in_edited_file
            #print "entire_list", entire_list
            missing=[x for x in entire_list if x not in list_of_samples_in_edited_file and [x] not in exceptions]
            #print "missing", missing
            print "missing from this file", sort_list_of_strings_that_are_actually_ints(missing)
            
    
    
#gets rid of all of the non-data columns at the top of the allele reports
def cleanup(data):
        rows_with_samples=[]
        #find the rows with a non empty column B and assemble a list with their indeces
        for row in data:
                #print "row", row
                if row and (row[1] or row[2]):
                        rows_with_samples.append(data.index(row))
        newdata=[]
        for a in range(0, len(rows_with_samples)):
                newdata.append(data[rows_with_samples[a]])
        return newdata

#takes a table (list of lists), converts empty values from the first row into their previous values
def fill_header(data):
        #print "data", data
        header=data[0]
        header.pop(0)
        header.pop(0)
        header.insert(0,'ID')
        header.insert(0,'Number')
        for y in range (2, len(header)-1):
                if header[y] and header[y+1]:
                        print "Might be a 1-column locus at ", header[y]
                        stop_hack=raw_input("You should probably stop and fix this") #figure out whether there's a more kosher way to kill the code here
        for x in range (0, len(header)): #header[0] and [1] should be blank
                if not header[x]:
                        header[x]=header[x-1]
        data.pop(0)
        data.insert(0,header)
        return data

#removes all columns that share labels with 3 or more other columns ---- and reports when one is removed that's on the list of ones that shouldn't be removed
def remove_bad_loci(data, keeplist):
        header=data[0]
        #print header
        counts=[]
#makes a list (counts) that counts the number of times an element appears in header
        for x in range (2, len(header)): #header[0] and [1] should be blank
                counts.append(header.count(header[x]))
        counts.insert(0,'')
        counts.insert(0,'')
        #print counts
#now make that list the first row of the data table
        data.insert(0,counts)
        #print data
#now transpose the whole table (easier to delete rows than columns)
        trans=transposed(data)
        #print trans
#now delete rows with value of 3 or greater as their first value (i.e. the ones where there are three or more alleles per individual)
        indeces_of_removals=[]
        for a in range(0, len(trans)):
                aslist=[]
                if trans[a][0]:
                        if trans[a][0] > 2:
                                aslist.append(trans[a][1]) #keeplist is a list of lists, so in order for the "in" comparison to work below, I needed to make trans[a][1] a list instead of a string
                                if aslist in keeplist:
                                        print trans[a][1] #this prints the locus name (of the loci that have more than three alleles, but should be nice and diploid)
                                indeces_of_removals.append(trans.index(trans[a]))

        while indeces_of_removals:
                trans.pop(indeces_of_removals[len(indeces_of_removals)-1])#-1?
                indeces_of_removals.pop()       

#now transpose what remains back, and return it
        trans2=transposed(trans)
        trans2.pop(0)
        return trans2

def get_unique(all_list): #gets the unique elements of a list and returns them as a list
        checked=[]
        for e in all_list:
                if e not in checked:
                        checked.append(e)
        return checked

def remove_nonscores(data, monolist):
        trans=transposed(copy.deepcopy(data))
        trans_copy=copy.deepcopy(trans)
        tracker=[]
        #print "trans", trans
        for a in monolist:
                #print "a", a
                for i, row in enumerate(trans):
                        #print "i at this point", i
                        if a[0] in row:
                                #print a[0], "is in row", i
                                tracker.append(i)
        #print "tracker", tracker
        tracker.sort()
        #print "sorted", tracker
        while len(tracker)>0:
                #print "tracker", tracker
                del trans_copy[tracker[len(tracker)-1]]
                del tracker[len(tracker)-1]
                               
                                
        new_data=transposed(trans_copy)
        #print "with nonscores removed", new_data
        return new_data

#takes a list of loci that were missing in each file (so, a locus is EXPECTED to be missing from every file except for the one it belongs in). This module takes that list, counts how many times each locus was missing, and returns a list of loci that were missing the maximal number of times
def most_missingest (cumulative, keeplist, pipeline_directory):
        unique_cumulative=get_unique(cumulative)
        #print "The loci that get scored", unique_cumulative
        missing=[]
        for w in range(0, len(keeplist)):
                if keeplist[w][0] not in unique_cumulative:
                        missing.append(keeplist[w])
        saveCsv(missing, 'missing.csv', pipeline_directory)


def rename_row (row, delimiter):#row must be a list of strings
        for x in range(1, len(row)):
                cut=row[x].find(delimiter)
                if cut > 0:
                        row[x]=row[x][:cut]
                else:#9_24_2011 this was added
                        cut=row[x].find('_')
                        if cut>0:
                                row[x]=row[x][:cut]#9_24_2011 this was added
                        
                        else:#if the above doesn't work, undo an indent for this else statement
                                print row[x], " couldn't be renamed-check manually. Please input a name you would like to call it (i.e. an ID number)"
                                new_name=raw_input()
                                row[x]=new_name
                        ##...
        return row


def index_of_longest_row(data): #can easily be modified to return the actual value instead of the index
        i_of_longest=0
        longest=0
        for x in range (0, len(data)):
                if len(data[x])>longest:
                        longest=len(data[x])
                        i_of_longest=x
        return i_of_longest

def index_of_row_with_x_in_slot_y(data, x, y): #data is a list of lists of strings
        #print x
        #print y
        for a in range(1, len(data)):
                #print data[a]
                #print data[a][y]
                if data[a][y]==x:
                        return a
        else:
                return "Your index_of_row_with_x_in_slot_y module had a problem"

def combine_allele_report_pipeline_dict(pipeline_directory, keeplist_address, rejected_samples_address):
        #The first allele report on the allele report list must contain the most comprehensive list of samples (i.e., if you have two DNA plates for the same panel of loci in two different GeneMarker files for one panel you scored the two plates together in another panel, use the GeneMarker file of the latter as the first in your list.
        files_to_combine=readCsv_input(pipeline_directory+'edited_file_list')

        address_to_directory=pipeline_directory
        keeplist=readCsv_input(keeplist_address)
        exceptions=readCsv_input(rejected_samples_address)
        #print "exceptions", exceptions
        whole_file=readCsv_input(address_to_directory+files_to_combine[0][0]) #I think the point of whole_file is to give a size to populate the list below
        #print "whole file", whole_file
        list_of_sample_names=[]
        for x in range(1, len(whole_file)):
                #print "whole_file[x][1]", whole_file[x][1]
                if [whole_file[x][1]] not in exceptions:
                        #print "this is not in exceptions, supposedly", whole_file[x][1]
                        list_of_sample_names.append(whole_file[x][1])
        #print "list of sample names", list_of_sample_names

        #populate a dicitonary called samples with keys that are the sample names with a list of locus names, each of which is itself a list that will later have genotypes appended to it
        samples={}
        location_of_keeplist=keeplist_address
        for x in range (0, len(list_of_sample_names)):
                #print "list_of_sample_names[x]", list_of_sample_names[x]
                samples[list_of_sample_names[x]]=readCsv_input(location_of_keeplist)
        #print "samples after filling", samples
        #print "the keys in samples", samples.keys()
        
        address_of_rejected_samples=rejected_samples_address
        exceptions=readCsv_input(address_of_rejected_samples)
        #print "exceptions", exceptions

        for w in range(0, len(files_to_combine)): #for each file in the list
            filename=address_to_directory
            filename+=files_to_combine[w][0]
            print 'file',filename
            temp=readCsv_input(filename) # read the file
            #print 'data from that file',temp
            #temp_trans=transposed(temp)
            target_index=0
            for z in range (1, len(temp)):
                #print "z", z
                for y in range (2, len(temp[0])):
                    #print "y", y
                    #print type(temp[z][1])
                    #print "temp[z][1] should be a sample name and a key in the samples dict",temp[z][1]
                    #print temp[z][1].find('B')
                    if (temp[z][1].find('blank')<0 and temp[z][1].find('Blank')<0 and [temp[z][1]] not in exceptions):#as long as the individual's name doesn't contain the word "blank", as in Blank or Blank_someNumber
                        #print "there should never be a blank here. Is there? :", temp[z][1] #the first "and" might need to change to an "or"
                        #print len(samples[int(temp[z][1])]) #might produce an error
                        for i in range (0, len(keeplist)):# I think is analogous #len(samples[int(temp[z][1])])
                                    #print "data [0][y]",temp[0][y]
                                    #print samples[int(temp[z][1])][i] #might produce an error
                                    #print "temp[0][",y,"]", temp[0][y]
                                    #print "temp[",z,"][1]", temp[z][1]
                                    #print "i", i
                                    #print "samples[temp[",z,"][1]][",i,"]", samples[(temp[z][1])][i]
                                    if temp[0][y] in samples[((temp[z][1]))][i]:# aka if the locus name is in this particular list within this particular sample #took out [i] #just added in the int here; might not be a good idea
                                        #print "this loop gets entered"
                                        #print int(temp[z][1])
                                        #print "i", i
                                        target_index=i
                                        #print "before",samples[(temp[z][1])][target_index]
                                        if (temp[z][1] not in exceptions  ): #this may not matter anymore because of the "exceptions" generated in the upstream script
                                            #print "target index", target_index
                                            #print "temp[z][y]", temp[z][y]
                                            if (temp[z][y]!='OL'and temp[z][y]!=''): # removed "and temp[z][y]!=''" from if statement
                                                if len(samples[(temp[z][1])][target_index])>3: #7_9_2012 I added the following script, which should collapse to diploid all samples that were tetraploid (likely because the samples were being repeated?) This, clearly, also deals with the problem of having two samples with the same names and different genotypes at a locus. It does so in such a way that missing genotypes are replaced with scored genotypes and homozygous genotypes are replaced with heterozygous ones.
                                                        #print "data [0][y]",temp[0][y]
                                                        #print samples[int(temp[z][1])][i] #might produce an error
                                                        #print "temp[0][",y,"]", temp[0][y]
                                                        #print "temp[",z,"][1]", temp[z][1]
                                                        #print "samples[temp[",z,"][1]][",i,"]", samples[(temp[z][1])][i]
                                                        old_genotype=samples[(temp[z][1])][target_index]                                                        
                                                        #print "old_genotype at this slot", old_genotype
                                                        new_genotype=[temp[z][y]]
                                                        #print "new_genotype at this slot", new_genotype
                                                        if '**' in old_genotype and '**' not in new_genotype:
                                                                samples[(temp[z][1])][target_index]=[ x for x in samples[(temp[z][1])][target_index] if "**" not in x ]                                                                
                                                                #print "after removal", samples[(temp[z][1])][target_index]
                                                        elif new_genotype[0] in old_genotype:
                                                                #print "had to remove one"
                                                                samples[(temp[z][1])][target_index].remove(new_genotype[0])
                                                                #print "after removal: ", samples[(temp[z][1])][target_index]
                                                        else:
                                                                #print "you haven't seen this problem yet", samples[(temp[z][1])][target_index]
                                                                Purge_Polyploid_MisScore_from_Locus(samples[(temp[z][1])][target_index])
                                                                #print "after ", samples[(temp[z][1])][target_index]
                                                samples[(temp[z][1])][target_index]=samples[(temp[z][1])][target_index]+[temp[z][y]]#may possibly need to deal with bracketing issues, depending on whether lists are appended to the rows, or strings are appended to the end of the single list in each row of each "keeplist"
                                                #print "after",samples[(temp[z][1])][target_index]
                                            else:
                                                #print "z",z,"y",y
                                                samples[(temp[z][1])][target_index]=samples[(temp[z][1])][target_index]+['**']
    
        #print "before re-purging", samples
        #testDict(samples)
        samples=Purge_Polyploid_MisScores(samples)
        DictOfLists_toTable(samples, pipeline_directory, rejected_samples_address)
        keys_as_numbers_sorted=GenePop_Output(samples, pipeline_directory, rejected_samples_address)
        return [list_of_sample_names, keys_as_numbers_sorted]

def any_line_that_starts_with_this(string, linelist):
        #print "string being passed", string
        any_trues=0
        for line in linelist:
                #print "line", line
                #print "type(line)", type(line)
                if line.startswith(string+','):
                        #print "this line is a hit:", line
                        any_trues+=1
        return any_trues

def sensitivity_analysis_HWE(pipeline_directory, sample_names_list, keys_as_numbers_sorted, set_demarcation_file_address):
        #this will take the complete GenePop format file with all sets of data and return GenePop formatted files with certain designated sets removed
        #print "sample_names_list_as_passed_to_sensivitiy_analysis", sample_names_list
        #print "keys_as_numbers_sorted_as_poassed_to_sensitivity_analysis", keys_as_numbers_sorted
        print "Beginning sensitivity analysis. This analysis visualizes the differences in deviation from Hardy Weinberg Equilibrium that the dataset displays with the removal of each plate as well as with each plate being used exclusively\n"
        
        linelist_input=pipeline_directory+'/final_output_GenePop.txt' 
        linelist=open(linelist_input).readlines()

        set_demarcation=readCsv_input(set_demarcation_file_address)
        set_demarcation_2=copy.deepcopy(set_demarcation)
        #print "set_demarcation_2", set_demarcation_2

        set_number=raw_input("How many DNA plates of samples are there? E.g., 7")
        #set_size=raw_input("How many samples per plate? E.g., 96")
        #print "sample_names_list", sample_names_list
        
        linelist = filter(lambda x: not x.isspace(), linelist)
        #print "linelist", linelist
        #print "linelist[0]", linelist[0]

        #get all of the loci and population information
        potential_POPs=['POP','Pop','pop']
        index_values=[]
        for x in range(0, len(potential_POPs)):
                #print "str(potential_POPs[x])",str(potential_POPs[x])
                #print "type", type(x)
                #print linelist.find(x) #will throw an error
                for y in range(0, len(linelist)):
                        if linelist[y].startswith(potential_POPs[x]):
                                index_values.append(y)
                  #try:
                 #       print "testing try"
                  #      print "linelist.find(potential_POPs[x])", linelist.find(potential_POPs[x])
                   #     index_values.append(linelist[linelist.find(potential_POPs[x])])
                #except:
                 #       print "exception happens"
                  #      pass
        #print "index_values", index_values
        top_part=linelist[0:min(index_values)+1]
        #print "top_part", top_part
        
        #create files excluding a particular set
        #print "excluding a particular set:"
        for x in range(0, int(set_number)): #make sure the end indexes right here
                flag=1
                set1_removed=open(pipeline_directory+'set'+str(x+1)+'_removed_GenePop_ready.txt', "w")
                #print "the file: ", pipeline_directory,'set',str(x+1),'_removed_GenePop_ready.txt'
                #for b in top_part:
                 #       set1_removed.write(b)
                #index=(x*int(set_size)) #deleted +1
                #back_index=((x+1)*int(set_size))-1 #just added -1
                index=set_demarcation[0][0]
                back_index=set_demarcation[1][0]
                #print "index", index
                #print "back_index", back_index
                
                for line in linelist:
                        #print "line", line
                        if line.startswith(index+','): #how to deal with getting the possibility that this sample is missing
                                #print index, " is the beginning of a line"
                                #print "line", line
                                flag=0
                        if line.startswith(back_index+','):
                                #print back_index, " is the beginning of a line"
                                #print "line", line
                                flag=1
                        if flag:#and not line.startswith("466")
                                set1_removed.write(line)
                set1_removed.close()
                del set_demarcation [0]
                del set_demarcation [0]
                #print "set_demarcation after pop", set_demarcation
                
        #create files consisting of only a particular set
        #print "only a particular set:"
        #print "top_part again:", top_part
        for x in range(0, int(set_number)): #make sure the end indexes right here
                flag=0
                set1_removed=open(pipeline_directory+'set'+str(x+1)+'_only_GenePop_ready.txt', "w")
                #print "the file: ", pipeline_directory,'set',str(x+1),'_only_GenePop_ready.txt'
                for b in top_part:
                        set1_removed.write(b)
                #index=(x*int(set_size))#removed +1
                #back_index=((x+1)*int(set_size))
                index=set_demarcation_2[0][0]
                back_index=set_demarcation_2[1][0]
                #print "index", index
                #good_starting_sample_for_the_set=sample_names_list[index]
                for line in linelist:
                    if line.startswith(index+','): #how to deal with getting the possibility that this sample is missing
                        #print index, " is the beginning of a line"
                        flag=1
                    if line.startswith(back_index+','):
                        #print back_index, " is the beginning of a line"
                        flag=0
                    if flag:#and not line.startswith("466")
                        set1_removed.write(line)
                set1_removed.close()
                del set_demarcation_2 [0]
                del set_demarcation_2 [0]
                #print "set_demarcation_2 after pop", set_demarcation_2

        
        print "Master files has been split and the individual sets and their respective removals from the master set are ready for HWE analysis"

def sort_list_of_strings_that_are_actually_ints(list_of_strings):
        #requires that they're actually ints. This will throw out the ones that aren't (but throws an error)
        #print "list_of_strings", list_of_strings
        new_list=[]
        for x in range (0, len(list_of_strings)):
                try:
                        new_list.append(int(list_of_strings[x]))
                except:
                        print "Did not include ", list_of_strings[x], " in the sorted list"
        return sorted(new_list)
                
def sortedDictValues(dictionary):
                keys=dictionary.keys()
                keys.sort()
                return map(dictionary.get, keys)
def genePop_to_rhh(pipeline_directory, rejected_samples_address):
        data=readCsv_input(pipeline_directory+'final_output_GenePop.txt')
        #print "data", data

        if (['POP'] in data):
                top=data.index(['POP'])+1
        else:
                if (['Pop'] in data):
                        top=data.index(['Pop'])+1
                else:
                        if (['pop'] in data):
                                top=data.index(['pop'])+1
                        
        #print "top", top
        data_trimmed=data[top:len(data)]
        #print "data_trimmed off top", data_trimmed
        #print data_trimmed[0]

        for x in range(0, len(data_trimmed)):
            if (len(data_trimmed[x])>1):
                #print type(data_trimmed[x][1])
                data_trimmed[x][1]=data_trimmed[x][1].split()
                #print "data[x] after splitting",data_trimmed[x][1]
                for z in range (0, len(data_trimmed[x][1])):
                    data_trimmed[x][1][z]=split(data_trimmed[x][1][z])
                    #print "data[",x,"][", 1,"]", "[",z,"]"
                #data_trimmed[x][1]=split(data_trimmed[x][y])

        #print "data_trimmed",data_trimmed
        for entry in data_trimmed:
                for genotype in entry[1]:
                        for i, allele in enumerate(genotype):
                                #print allele
                                if allele in ['00','000']:
                                        genotype[i]='0'
                                        #print genotype 

        #print "data_trimmed", data_trimmed
        #print type(data_trimmed)
        
        spaces_removed=remove_spaces_from_list_recursive(data_trimmed)
        #print "spaces_removed", spaces_removed
        rhh_output(data_trimmed, pipeline_directory)

        print "genePop data has been converted to rhh format and is ready for the rhh script"

def remove_spaces_from_list_recursive(the_list): #can deal with lists of lists of lists etc...
        if isinstance(the_list, basestring):
            #print "base case happens", the_list
            return the_list.replace(" ","")
        else:
            #print "base case didn't happen for", the_list
            return map(remove_spaces_from_list_recursive, the_list)                       
                                        

def convert_to_g2 (pipeline_directory, keeplist_address):
        #this script will take a file in GenePop format and will convert it into a file that is ready to be submitted as RMES input (with the possible exception of having to deal with mac vs pc issues, and perhaps a space at the end)

        #Assumptions: each row is the same length

        #get the GenePop text file, and convert it to csv format
        address=pipeline_directory+'final_output_GenePop.txt'
        txt_to_csv(address)
        data=readCsv_input(address+'_csv_converted.csv')
        #print "data", data


        #extract the important parameters of the GenePop file. Attempting to generalize to GenePop files with more than one population
        number_of_populations=data.count(['POP'])+data.count(['Pop'])+data.count(['pop'])+data.count([['POP']])+data.count([['Pop']])+data.count([['pop']])
        print "number of populations", number_of_populations

        #now get a list of population names. This is only complicated because there could be more than one, and it's hard to distinguish them in a list
        #assumes that the GenePop file has more than one population name if it has more than one population. This seems not to be the case, but RMES suggests that there might. I'll create a way to manually add more population names (may not be useful). There is potential that this will just create a list of the same population name over and over again. If there is an error here and it is actually a problem for your dataset, please email mark.aaron.fisher@gmail.com
        population_names=[]
        keeplist=readCsv_input(keeplist_address)
        #print "keeplist", keeplist
        #When you made the GenePop conversion module, you changed some the locus names to avoid having spaces in them. You replaced them with underscores. Here, you'll put the spaces back where the underscores are. *This is why your locus names can't have underscores already!*
        keeplist_without_spaces=keeplist
        for x in range(0, len(keeplist_without_spaces)):
            #print "before replace",keeplist_without_spaces[x][0]
            keeplist_without_spaces[x][0]=keeplist_without_spaces[x][0].replace(' ', '_')
            #print "after replace", keeplist_without_spaces[x][0]
        #print "keeplist_without_spaces", keeplist_without_spaces
        #still getting population names...
        name_counter=0
        for x in range (0,number_of_populations):
            #print "x", x
            for y in range (0, len(data)):
                data=readCsv_input(address+'_csv_converted.csv')
                #print data[y]
                #print type(data[y])
                if (len(data[y])>0):                        
                        if (([data[y][0]] not in keeplist) and ([data[y][0]] not in keeplist_without_spaces) and (data[y][0] not in ('POP', 'Pop', 'pop', ['POP'], ['Pop'], ['pop'])) and (name_counter<x+1)): #the name_counter here tracks whether the potential population name occurs before the word 'POP' in this particular instance of finding a population name. In other words, it ensures that population_name.append() only happens once per loop. This avoids lines with genotypes or locus names being reported as population names
                                population_names.append(data[y][0])
                                name_counter+=1
        print "These are what I have as your population names:", population_names
        done=raw_input("Is that all of them? y/n")
        while done != ('y' or 'yes'):
            new_pop_name=raw_input("Please enter the next population name.")
            population_names.append(new_pop_name)
            done=raw_input("Is that all of them? y/n")
        print "population_names", population_names


        #now, start to split the genotype strings and convert them to 1 for het, 0 for homo, -99 for missing score
        for y in range (0, len(data)):
            #print "y in split module", y
            if (len(data[y])>0):
                if ([data[y][0]] not in keeplist and ([data[y][0]] not in keeplist_without_spaces) and (data[y][0] not in population_names) and (data[y][0] not in ('POP', 'Pop', 'pop', ['POP'], ['Pop'], ['pop']))):
                    #print "data[",y,"][0]",data[y][0]
                    data[y][0]=data[y][0].split(' ')
                    #print "before", data[y][0]
                    for z in range (1, len(data[y][0])):
                        data[y][0][z]=assign_zygosity(data[y][0][z])
                    #print "after zygosity assignment", data[y][0]
            

        #print "data at this point", data

        #create a csv file that contains ids, then each locus and the genotype at that locus for each sample as either 1, 0, or NA
        hets_and_homos_cleanup(copy.deepcopy(data), population_names, pipeline_directory) #this takes data at this point and converts it into a csv file that can later be used for calculate_MLH. It's a completely internal module.

        #create a csv file that contains ids in one column, and MLH in the other
        calculate_MLH(readCsv_input(pipeline_directory+'heterozygotes_and_homozygotes.csv'), pipeline_directory)

        
        #get the number of loci and the number of samples in each population
        number_loci=[]
        sample_number_counter=[]

        for x in range (0, number_of_populations):
            POP_counter=0
            #print "x",x    
            for y in range (0, len(data)):
                #print "y", y
                #print "POP_counter at this point", POP_counter
                if (POP_counter==0 and ([data[y][0]] in keeplist or [data[y][0]] in keeplist_without_spaces)):
                    if (len(number_loci)<x+1):
                        number_loci.append(1)
                        #print "made a new locus element"
                    else:
                        number_loci[x]=number_loci[x]+1
                        #print "added a new locus to the current population"
                if (POP_counter==x+1 and ([data[y][0]] not in keeplist and [data[y][0]] not in keeplist_without_spaces) and (data[y][0] not in ('POP', 'Pop', 'pop', ['POP'], ['Pop'], ['pop'])) and (data[y][0] not in population_names)):
                     if (len(sample_number_counter)<x+1):
                        sample_number_counter.append(1)
                        #print "made a new sample number element"
                     else:
                         #print sample_number_counter[x]
                         sample_number_counter[x]=sample_number_counter[x]+1
                         #print "added a new sample to the current population"
                if (data[y][0] in ('POP', 'Pop', 'pop', ['POP'], ['Pop'], ['pop'])):
                    #print "POP_counter gets added to"
                    POP_counter+=1
        
        #now create the output file
        #print "data", data       
        output_name=pipeline_directory+'RMES_ready.txt'
        saveFile=open(output_name,'w')
        saveFile.write(str(number_of_populations))
        saveFile.write('\n')
        for x in range (0, number_of_populations):
            #print "x",x
            saveFile.write(population_names[x])
            saveFile.write('\n')
            #print sample_number_counter[x]
            #print "type that above is", type(sample_number_counter[x])
            saveFile.write(str(sample_number_counter[x]))
            saveFile.write('\n')
            saveFile.write(str(number_loci[x]))
            saveFile.write('\n')
        for y in range (0, len(data)):
            #print "y",y
            if (data[y][0] not in keeplist and data[y][0] not in keeplist_without_spaces and data[y][0] not in population_names and data[y][0] not in ('POP', 'Pop', 'pop', ['POP'], ['Pop'], ['pop']) and [data[y][0]] not in keeplist and [data[y][0]] not in keeplist_without_spaces and [data[y][0]] not in population_names and [data[y][0]] not in ('POP', 'Pop', 'pop', ['POP'], ['Pop'], ['pop'])):
                #print "data[",y,"][0]", data[y][0]
                #print "type", type(data[y][0])
                if len(data[y][0])>1:
                    for z in range (1, len(data[y][0])):
                        #print "data[",y,"][0][",z,"]", data[y][0][z]
                        saveFile.write(str(data[y][0][z])+' ')
                    saveFile.write('\n')                               
                else:
                    #print "this is your list that should be a string:", data[y][0][0]
                    saveFile.write(data[y][0][0])
                    saveFile.write('\n')
        saveFile.close()
        print "Your GenePop file has been converted to an RMES-ready file called, ",output_name

def PypeRify(r_as_text_address): #gets a little buggy with pounded comments on the same lines as R code
        warning=raw_input("Just remember that your R code should have no pounded comments before entering this module, or else it won't work")
        linelist=open(r_as_text_address).readlines()
        saveFile=open(raw_input("Address and name-including '.txt'-of the file you want to save"), 'w')
        for line in linelist:
                #print "type(line)", type(line)
                #print "line before", line
                if len(line)>0:
                        line=line.replace('"','\'')
                        line=insert(line,'print r("',0)
                        line=insert(line,'")', len(line)-1)
                        #print "line after", line
                        saveFile.write(line)
                        #print line
                        #print "type", type(line)
                        #def insert(original, new, pos)

        saveFile.close()

def Purge_Polyploid_MisScores(dictOfLists):
    for key,lst in dictOfLists.items():
        for locus in lst:
            while len(locus)>3:
                #print "locus %s" %locus
                if '**' in locus:
                    locus.remove('**')
                else:
                    for allele in locus:
                            #print "allele %s" %allele
                            if locus.count(allele)>2:
                                    locus.remove(allele)
                            elif locus.count(allele)>1 and len(locus)>3:
                                locus.remove(allele)
                #print "still in the loop"
                
    return dictOfLists

def Purge_Polyploid_MisScore_from_Locus(Locus_list):
        if len(Locus_list)>3:
                try:
                        Locus_list=Locus_list.remove('**')
                except:
                        for z in Locus_list:
                                if Locus_list.count(z)>2:
                                        #print "removed ",z," at Locus_list", Locus_list
                                        Locus_list.remove(z)
                                        #print "now, it looks like", Locus_list
                if len(Locus_list)>3:
                        print "The Length is still greater than 3! at Locus_list", Locus_list
        return Locus_list
                                        

        
def Incorporate_HL_from_rhh(pipeline_directory, address_of_acceptor_file):
        #takes the rhh MLH output (ID, IR, SH, and HL columns), and converts it into a csv file on the desktop called csv_converted.csv
        txt_to_csv(pipeline_directory+"Rhh_test_output.txt")


        #makes said csv file the donor file
        #print "the prompt you're about to get is for the donor file-which is the file that contains the ID, IR, SH, and HL columns from the rhh script. It might be this file: /Users/markfisher/Desktop/pipeline/Rhh_test_output.txt_csv_converted.csv"
        donor=readCsv_input(pipeline_directory+'Rhh_test_output.txt_csv_converted.csv')
        #print "donor",donor, "done now"
        acceptor=readCsv_input(address_of_acceptor_file)
        #print "acceptor", acceptor

        acceptorTrans=transposed(acceptor)
        #print "acceptorTrans", acceptorTrans

        #print donor
        tally=0
        for i in range (1,len(donor)):
            if donor[i][0] in acceptorTrans[0]:
                tally+=1
                acceptor[acceptorTrans[0].index(donor[i][0])].append(donor[i][1])
                acceptor[acceptorTrans[0].index(donor[i][0])].append(donor[i][2])
                acceptor[acceptorTrans[0].index(donor[i][0])].append(donor[i][3])
        #print acceptor
        #print "tally", tally

        #add IR, SH, HL columns
        acceptor[0].append('IR')
        acceptor[0].append('SH')
        acceptor[0].append('HL')
        #print "acceptor before", acceptor
        acceptor=fill_blanks_of_data_table(acceptor)
        acceptor=make_sure_NA_are_numeric_compliant(acceptor)
        acceptor=remove_spaces_from_header(acceptor)
        saveCsv(acceptor,"MLH_output.csv", pipeline_directory)

        print "IR, HL, and SH values have been incorporated into the fitness data file (called MLH_output.csv)"

def remove_spaces_from_header(list_of_lists):
        for entry in list_of_lists[0]: #assumes the headers are in the top row of list_of_lists
                entry.replace(" ","")
                entry.replace("-","_")
        return list_of_lists

def fill_blanks_of_data_table(data):
        #print data
        #assumes data is a list of lists
        max_len=[]
        for x in data:
                max_len.append(len(x))
        the_max=max(max_len)
        #print "the max is", the_max
        for x in data:
                #print len(x)
                if len(x)<the_max:
                        #print "this happens"
                        while len(x)<the_max:
                                x.append('NA')
                for i,y in enumerate(x):
                        #print "y", y
                        if y=='' or len(y)==0:
                                #print "this happens", y
                                x[i]='NA'
        #print "data after filling blanks", data
        return data

def clean_and_save(filename, index, filenames, destination):
    flag=0
    #print "filename inside clean_and_save module", filename
    output_file=open(destination+filenames[index]+'_output.txt','w')
    linelist=open(filename).readlines()
    #print "linelist", linelist
    for line in linelist:
        if line.startswith('locus'):
            flag+=1
        if line.startswith(' Chi2'):
            flag+=1
        if flag==1 and not line.startswith('-') and not line.startswith(' Chi2') and not line.startswith('locus') and not line.startswith('All'): # and not line.startswith('')
            #print "this happens"
            #print line
            output_file.write(line)
    output_file.close()
    name_to_be_split=destination+filenames[index]+'_output.txt'
    #print "name_to_be_split", name_to_be_split
    txt_to_csv_space_split_version(name_to_be_split)
        #print "flag", flag
    ##remove_two_word_locus_names(name_to_be_split+'_csv_converted.csv') #I hope that getting rid of this was ok!


def GenePop_HWE_to_Pvalue_csv(address_of_sensitivity_GenePop_directory):
        #if the directory doesn't exist, make it
        try:
                os.mkdir(address_of_sensitivity_GenePop_directory)
        except:
                pass

        #make a subdirectory for the output
        try:
                os.mkdir(address_of_sensitivity_GenePop_directory+'Converted_to_csv') # raises an error if the directory already exists
        except:
                pass
        
        #wait for the user to put the GenePop html text files into the relevant directory
        answer=raw_input("Have you run GenePop webversion on your GenePop files (just the master one, or all of the subfiles for sensitivity analysis) and saved the output as text files into "+address_of_sensitivity_GenePop_directory+"? y/n")
        if (answer=="y" or answer=="yes"):

                #get a list of all of the .txt files in it
                #print "address_of_sensitivity_GenePop_directory: ",address_of_sensitivity_GenePop_directory
                file_names=os.listdir(address_of_sensitivity_GenePop_directory)
                #print "file names list: ", file_names
                #filter the list so that it's only txt files
                file_names_filtered=[x for x in file_names if '.txt' in x]
                #print "file_names_filtered", file_names_filtered
                
                for y in range(0, len(file_names_filtered)):
                    #print "this happens 2"
                    clean_and_save(address_of_sensitivity_GenePop_directory+file_names_filtered[y],y,file_names_filtered, address_of_sensitivity_GenePop_directory+'Converted_to_csv/')
        return address_of_sensitivity_GenePop_directory+'Converted_to_csv/'
def pvalues_of_all(sensitivity_csv_subdirectory_address):
        file_names=os.listdir(sensitivity_csv_subdirectory_address)
        file_names_filtered=[x for x in file_names if '.csv' in x]
        data=[['Loci']]
        #print "file_names_filtered", file_names_filtered

        #make a pvalues all directory (otherwise, Pvalues_of_all itself will be counted among the csv files in this directory)
        try:
                os.mkdir(sensitivity_csv_subdirectory_address+'Pvalues_combined')
        except:
                pass
        #pvalues_file=open(sensitivity_csv_subdirectory_address+'Pvalues_of_all.csv','w')

        #create header row
        file_names_as_headers=[x[0:x.find('.')] for x in file_names_filtered]
        #print "file_names_as_headers", file_names_as_headers
        for x in file_names_as_headers:
                data[0].append(x)
        #print "data", data

        #create locus name column
        locus_getter=readCsv_input(sensitivity_csv_subdirectory_address+file_names_filtered[0]) #all of the files should have the same locus names
        for a in locus_getter: #had to get rid of an empty entry at the end of the file
                #print "a", a
                if not a:
                        #print "no a here", a
                        locus_getter.pop(locus_getter.index(a))
        #print "locus_getter", locus_getter
        locus_names=transposed(locus_getter)[0]
        #print "locus_names", locus_names
        for x in locus_names:
                data.append([x])
        #print "data", data

        #add values from the second column of each csv file into data
        for x in file_names_filtered:
              current_file=readCsv_input(sensitivity_csv_subdirectory_address+x)
              #print "current_file", current_file
              for y in range(0, len(current_file)-1):#unclear to me why this should be -1, especially since I fixed the empty entry problem in locus_getter (see above)!
                      #print "y", y
                      #print "data[y+1]", data[y+1]
                      #print "current_file[y][1]", current_file[y][1]
                      data[y+1].append(current_file[y][1])
        #print "data", data
        #print "sensitivity_csv_subdirectory_address", sensitivity_csv_subdirectory_address
        saveCsv(data, 'Pvalues_of_all.csv', sensitivity_csv_subdirectory_address+'Pvalues_combined/')       
        

def get_bifurcation(data, pipeline_directory, rejected_samples_address):
        #print "data in bifurcation function", data
        data_as_dict=Table_to_DictOfLists_General(data) # remember that this module assumes that the first column is IDs (and thus unique)
        exceptions=readCsv_input(rejected_samples_address)
        #print "data_as_dict", data_as_dict
        try:
                os.mkdir(pipeline_directory+'Data_subsets')
        except:
                pass
        
        #print all of the column options
        print "Your trait column headers are as follows (BE MINDFUL OF SPACES!:\n", data[0]
        print "Your trait values include such values as (BE MINDFUL OF SPACES!:\n", data[1]
        #get an array from the user of column
        new_value='nope'
        print "Which traits and which values of these traits do you want to exclude from this analysis? Please enter as two values separated by a comma (NO SPACES!). Example: V1,0. Simply type 'done' when done."
        traits_and_values=[]
        while new_value not in ['done','Done']:
                #print "new_value", new_value, "type", type(new_value)
                new_value=raw_input()
                split_new_value=new_value.split(',')
                split_new_value_remove_spaces=[x.replace(' ','') for x in split_new_value] #can't remember why I originally did this, but it creates a problem when the trait mismatches the trait list, which contains spaces. Removed for now 7.17.2012
                print "temp split_new_value_remove_spaces[0]", split_new_value_remove_spaces[0]
                the_trait=split_new_value_remove_spaces[0]
                for i,row in enumerate(transposed(data)):
                        if row[0]==the_trait:
                                the_row=i
                                #print the_trait, "=", row[0]
                                
                #print "type", type(split_new_value_remove_spaces[0])
                #print "data", data
                #print "data.index(split_new_value_remove_spaces[0])", data.index(split_new_value_remove_spaces[0])
                if (split_new_value_remove_spaces[0] not in data[0] and split_new_value_remove_spaces[0] not in ['done','Done']):
                        print "You have entered a trait that is not found in the trait-containing spreadsheet. Please re-enter, being mindful of spaces and typos."
                elif (len(split_new_value_remove_spaces)!=2 and split_new_value_remove_spaces[0] not in ['done','Done']):
                        print "You may have forgotten to enter a trait value along with your trait"
                elif (len(split_new_value_remove_spaces)==2):
                        if (split_new_value_remove_spaces[1] not in transposed(data)[the_row] and split_new_value_remove_spaces[1] not in ['done','Done']):
                                        print "You have entered a trait value that is not found in the trait-containing spreadsheet. Please re-enter, being mindful of spaces and typos."
                        else:
                                traits_and_values.append(split_new_value) # changed from traits_and_values.append(split_new_value_remove_spaces)
                #print "traits_and_values", traits_and_values
        #remove the 'done' from traits_and_values
        #I don't think you have to do this anymore
        #traits_and_values=traits_and_values[0:len(traits_and_values)-1]
        print"traits_and_values final", traits_and_values

        data_copy=copy.deepcopy(data_as_dict)

        #remove header
        #print "before header removal", data_copy
        data_copy=remove_header_dict(data_copy)
        #print "after header removal", data_copy
        
        #print "data_as_dict", data_as_dict
        #for x in data_as_dict:
#                print "x all", x
        for x in data_as_dict:
                #print "x", x
                for y in data_as_dict[x]:
                        #print "y", y
                        for z in traits_and_values:
                                #print "z", z
                                #print "z[0]", z[0]
                                #print "y[0]", y[0]
                                #print "y[1]", y[1]
                                #print "y", y
                                if z[0]==y[0]:
                                        #print "this happens" #z[0]==y[1] or or  removed from if statement below
                                        if z[1]==y[1] or  (' '+z[1])==y[1] or [x] in exceptions: #will the z[0]==y[1] ever cause problems? Using this currently to eliminate the header row that keeps ending up in the output
                                                #print "I will remove this sample from the dictionary"
                                                try:
                                                        #print x, " is getting deleted"
                                                        del data_copy[x]
                                                        
                                                except:
                                                        print "Hopefully the script has already deleted ", x
                                                #print x
                                                #print z
                                                #print y
                                                
                                        
        #print "data_copy", data_copy
        try:
                os.mkdir(pipeline_directory+'Data_subsets/'+list_as_string([list_as_string(x) for x in traits_and_values])+'_removed_subset/')
        except:
                pass
        saveCsv(DictOfLists_toTable_General(data_copy, pipeline_directory), list_as_string([list_as_string(x) for x in traits_and_values])+'_removed_subset.csv', pipeline_directory+'Data_subsets/'+list_as_string([list_as_string(x) for x in traits_and_values])+'_removed_subset/')
        
        return [DictOfLists_toTable_General(data_copy, pipeline_directory), pipeline_directory+'Data_subsets/'+list_as_string([list_as_string(x) for x in traits_and_values])+'_removed_subset/']

def remove_header_dict(data):
    data_copy=copy.deepcopy(data)
    for x in data:
        count=0
        for y in data[x]:
            if y[0]==y[1]:
                count+=1
        if count==len(data[x]):
            del data_copy[x]
    return data_copy

def perform_MAM(data, pipeline_directory, rejected_samples_address):#assumes that exceptions have been removed at this point
        print "This portion of the pipeline evaluates minimum adequate models, to determine whether heterozygosity is an important predictor of a target fitness-associated trait.\n"
        #data_as_dict=Table_to_DictOfLists_General(data)
        #exceptions=readCsv_input(rejected_samples_address)
        
        #print "data_as_dict", data_as_dict
        try:
                os.mkdir(pipeline_directory+'Regressions')
        except:
                pass

        saveCsv(data,'temp.csv',pipeline_directory+'Regressions/')
        #print "Make sure the temp.csv file gets completely overwritten each time"

        #print all of the column options
        print "Your trait column headers are as follows (BE MINDFUL OF SPACES!:\n", data[0]
 

        #get the relevent link function for the response variable (and predictor variables?)
        family=raw_input("What sort of error distribution and link function do you expect to use in your models? Your choices are: binomial, gaussian, Gamma, inverse.gaussian, poisson, quasi, quasibinomial, and quasipoisson. Note that gaussian denotes normal linear regression. Please type one of them and press enter")
        while family not in ['binomial', 'gaussian', 'Gamma', 'inverse.gaussian', 'poisson', 'quasi', 'quasibinomial', 'quasipoisson']:
                family=raw_input("You typed a link function that is not one of the options. Please type either binomial, gaussian, Gamma, inverse.gaussian, poisson, quasi, quasibinomial, or quasipoisson. Note that gaussian denotes normal linear regression.")

       
        #get the response variable from the user
        new_value='not done'
        print "What trait do you want as your response variable?"
        new_value=raw_input()
        while new_value not in data[0]:
                print "You have entered a trait that is not found in the trait-containing spreadsheet. Please re-enter, being mindful of spaces and typos."
                new_value=raw_input()
        response_variable=new_value

        #get the predictor variable(s) from the user
        predictor_variables=[]
        print "What trait(s) do you want as predictor variable(s)? Please hit RETURN between each entry. Simply type 'done' when done."
        new_value='not done'
        while new_value!='done' and new_value!='Done':
                #print "new_value", new_value, "type", type(new_value)
                new_value=raw_input()
                if (new_value not in data[0] and new_value not in ['done','Done']):
                        print "You have entered a trait that is not found in the trait-containing spreadsheet. Please re-enter, being mindful of spaces and typos."
                else:
                        predictor_variables.append(new_value) # changed from traits_and_values.append(split_new_value_remove_spaces)
                #print "predictor_variables", predictor_variables

        #remove the 'done' from predictor_varaibles
        predictor_variables=predictor_variables[0:len(predictor_variables)-1]
        print"predictor_variables", predictor_variables


        
        #the bulk of the work (needs to be new PypeR stuff)
        predictors_as_string=list_as_string(predictor_variables)
        #from pyper import *
        r=R(use_numpy=True)
        r.pipeline_dir=pipeline_directory
        r.response=response_variable
        r.predictors=predictor_variables
        r.predictor_filename=predictors_as_string
        r.fam=family
        print r("setwd(paste(pipeline_dir,'Regressions/', sep=''))")
        print r("data<-read.csv(paste(pipeline_dir,'Regressions/temp.csv', sep=''),header=T)")
        print r("data[data=='-']=NA")
        print r("data[data=='']=NA")
        #print r("data[data==NULL]=NA")
        ####THIS IS WHERE YOU SHOULD MAKE YOUR EXCEPTION SUBSET
        r.rejected_path=rejected_samples_address
        print r("rejected<-read.csv(rejected_path, header=F)")
        #print r("print(rejected)")
        print r("data<-subset(data, !(ID %in% rejected[,1]))") 
        #print r("print(data[c(predictors, response)])")
        print r("no.na.data<-na.omit(data[c(predictors,response)])")
        #print r("str(data)")
        print r("model <- glm(formula=as.formula(paste(paste(response,'~', sep=''),paste(predictors,collapse='+'), sep='')),family=fam, no.na.data)")
        CMDs1=["sink(file=paste(paste(pipeline_dir,'Regressions/', sep=''),paste(paste(paste(paste(fam,response,sep=''), 'on', sep=''), predictor_filename, sep=''),'_fullModel.txt', sep='')))", "print(summary(model))", "sink()"]
        CMDs2=["sink(file=paste(paste(pipeline_dir,'Regressions/', sep=''),paste(paste(paste(paste(fam,response,sep=''), 'on', sep=''), predictor_filename, sep=''),'.txt', sep='')))", "step(model)", "sink()"]
        CMDs3=["tiff(filename=paste(paste(pipeline_dir,'Regressions/', sep=''),paste(paste(paste(paste(fam,response,sep=''), 'on', sep=''), predictor_filename, sep=''),'_fullModel_assumptions_check.tiff', sep='')))", "layout(matrix(c(1,2,3,4),2,2))", "plot(model)", "dev.off()"]
        print r(CMDs1)
        print r(CMDs2)
        print r(CMDs3)

def list_as_string(list_of_strings):
        final_string=''
        for s in list_of_strings:
                final_string+=s
                final_string+='_'
        return final_string[0:len(final_string)-1]

def Table_to_DictOfLists_General(table): #assumes each header row with ID in first col and each other cell in top row populated with locus name (ie. doubled), missing data contains **, allele values are at least two digits
        samples={}
        table_trans=transposed(table)
        #print table_trans[0]
        #print len(table_trans[0])
        #print "table", table
        headerList=[[x] for x in table[0]]
        #print"headerList",  headerList
        for x in table_trans[0]:
                samples[x]=copy.deepcopy(headerList)

        
        for y in samples:#for each individual
                #go through the traits:
                for z in range(1, len(samples[samples.keys()[1]])+1):
                        #print "does len change?", len(samples[samples.keys()[1]])+1
                        #print "y", y
                        #print "z-1", z-1
                        #print "samples[y][z-1]",samples[y][z-1]
                        #print "table_trans[0].index(str(y))", table_trans[0].index(str(y))
                        #print "table[table_trans[0].index(str(y))]", table[table_trans[0].index(str(y))]
                        #print "table[0].index(samples[y][z-1][0])]", table[0].index(samples[y][z-1][0])
                        #print "table[table_trans[0].index(str(y))][table[0].index(samples[y][z-1][0])]", table[table_trans[0].index(str(y))][table[0].index(samples[y][z-1][0])]
                        samples[y][z-1].append(table[table_trans[0].index(str(y))][table[0].index(samples[y][z-1][0])])
        return samples 

def DictOfLists_toTable_General(dictOfLists, pipeline_directory):                
        #address_of_rejected_samples=rejected_samples_address
        data=[['Key']]
        #populate the top row
        for i in range (0, len(dictOfLists[dictOfLists.keys()[0]])):
                data[0].append(dictOfLists[dictOfLists.keys()[0]][i][0])
                
        #go through all of the other rows
        for x in range (0, len(dictOfLists)):
                #print "x", x
                data.append([])
                data[len(data)-1].append(dictOfLists.keys()[x])# make sure these samples match their genotypes
                for y in range(0,len(dictOfLists[dictOfLists.keys()[0]])):
                        #print "y", y
                        #if (len(dictOfLists[dictOfLists.keys()[x]][y])==3):
                                #print "this happens at x", x, " and y", y
                        data[len(data)-1].append(dictOfLists[dictOfLists.keys()[x]][y][1])#just changed this from x+1 to x for indexing errors after generalizing. 6.24.2012
                        #data[len(data)-1].append(dictOfLists[dictOfLists.keys()[x]][y][2])#see lijne above
        #saveCsv(data, 'final_output', pipeline_directory)
        return data

def Correlation_heat_map(pipeline_directory, rejected_samples_address): #data is a dataframe R object
        try:
                os.mkdir(pipeline_directory+'Correlations')
        except:
                pass
        r=R(use_numpy=True)
        print r("library(lattice)")
        r.pipeline_dir=pipeline_directory
        r.address_of_acceptor_file_r=pipeline_directory+'MLH_output.csv'
        print r("data_for_cor<-read.csv(address_of_acceptor_file_r)")
        r.rejected_path=rejected_samples_address
        print r("rejected<-read.csv(rejected_path, header=F)")
        print r("data_for_cor<-subset(data_for_cor, !(ID %in% rejected[,1]))") 

        
        #print r("print(data_for_cor)")
        print r("just_nums<-sapply(data_for_cor,is.numeric)")
        #print r("print(just_nums)")
        print r("just_nums_data<-data_for_cor[,just_nums]")
        print r("cor2<-cor(just_nums_data, use='pairwise')")
        #print r("print(cor2)")
        print r("rgb.palette <- colorRampPalette(c('blue', 'yellow', 'red'), space = 'rgb')")
        CMDs1=["pdf(file=paste(paste(pipeline_dir,'Correlations/', sep=''), '_heat_map.pdf',sep=''))", "levelplot(cor2, aspect='iso', scales=list(x=list(rot=90)), main='Correlation Of Fitness-Associated Triats', col.regions=rgb.palette(220), cuts=200, at=seq(-1,1,0.01))"]
        print r(CMDs1)
        print r("dev.off()")

def calculate_MLH(data, pipeline_directory):
        data_copy=copy.deepcopy(data)
        del data_copy[0]
        #print data_copy
        output=[]
        global_loci_typed_stats=[]
        for row in data_copy:
                #print "row", row
                row_copy=copy.deepcopy(row)
                del row_copy[0]
                het_count=0
                total_typed_count=0
                for genotype in row_copy:
                        #print "genotype", genotype
                        if genotype in ['1', '0']:
                                #print "a score happens"
                                total_typed_count+=1
                        if genotype == '1':
                                #print "a het happens"
                                het_count+=1
                #print "het_count", het_count
                #print "total_typed_count", total_typed_count
                if total_typed_count==0:
                        MLH='NA'
                else:
                        MLH=float(het_count)/float(total_typed_count)
                #print "MLH", MLH
                
                output.append([row[0],MLH])
                global_loci_typed_stats.append([row[0],total_typed_count])
        output.insert(0,[data[0][0],'MLH'])
        global_loci_typed_stats.insert(0,[data[0][0],'Loci_typed_per_individual'])
        #print output
        saveCsv(output,'MLH_calculated.csv', pipeline_directory)
        saveCsv(global_loci_typed_stats,'Loci_typed_per_individual.csv', pipeline_directory)


def incorporate_donor_into_acceptor(acceptor_address,donor_address,output_address):
    #output_address should be a text file        
    #this script is for incorporating genotype information for only samples that are suitable


    good_sample_list=readCsv_input(acceptor_address)

    #print good_sample_list
    data=readCsv_input(donor_address)

    #print data
    #data_trans=transposed(data)
    #print data_trans
    replacementCount=0
    #temp=[]
    #print "len(data_trans)=", len(data_trans)
    for x in range(0, len(good_sample_list)):
        #print "x", x
        for y in range(0,len(data)):#not exactly sure why this should be -1
            #print "y", y
            #rint "that",good_sample_list[x][0]
            #print "this",data[y][0]
            if(good_sample_list[x][0]==data[y][0]):
                replacementCount+=1
                #print data_trans[y]
                good_sample_list[x][0]=data[y]
        #        temp=data_trans.pop(y) #delete that row
                #print temp
       # for z in range(0,len(data_trans)-1):#not exactly sure why this should be -1
            #print "y", y
            #print data_trans[y]
        #    if(unhardy[x][0]==data_trans[z][0]):
         #       replacementCount+=1
                #print data_trans[z]
          #      temp=data_trans.pop(z) #delete that row
                #print temp

    #data=transposed(data_trans)
    saveFile= open(output_address, 'w')
    for w in range(0, len(good_sample_list)):
        saveFile.write(str(good_sample_list[w]))
        saveFile.write('\n')
    saveFile.close()
    #print good_sample_list
    #saveCsv(good_sample_list, 'genotypes_incorporated')
    #print replacementCount
    #print 'unhardy', unhardy
    #print "data_trans", data

def incorporate_donor_into_acceptor_csv (acceptor_address, donor_address, output_address, output_filename, add_donors_missing_in_acceptor):
        acceptor=readCsv_input(acceptor_address)
        #print "acceptor", acceptor
        donor=readCsv_input(donor_address)
        #print "donor", donor

        for i, row in enumerate(acceptor):
                #print "row", row
                for r in donor:
                        #print "r", r
                        if row[0] == r[0]:
                                for item in r[1:len(r)]:
                                        #print "item", item
                                        acceptor[i].append(item)
                if len(acceptor[i])<len(donor[0]):
                        for x in range (0,len(donor[0])-1):
                                acceptor[i].append('NA')

        if add_donors_missing_in_acceptor in ['y','yes','Yes','Y']:
                ##check if there's any unique stuff in donor and add that to the end of the new acceptor
                #print "acceptor", acceptor
                trans=transposed(acceptor)
                #print "trans", trans[0]
                for r in donor:
                        if r[0] not in trans[0]:
                                acceptor.append(r)

        #print "acceptor", acceptor

        ##save
        saveCsv(acceptor,output_filename,output_address)
                
def hets_and_homos_cleanup(data, population_names, pipeline_directory): #right now might only work for single populations
        #assumed data is a list of lists, the first element containing a population name, the next group containing locus names (each its own element), then 'POP' or something like this, then sample names with 1s, 0s, and -99s. These genotype-containing elements are the only ones that contain more than one element themselves

        #split data into two lists of lists, the first with population and locus names and pop

        non_sample_info=[]
        sample_info=[]
        for row in data:
                #print "row", row
                #=print "type(row)", type(row[0])
                #so, the sample_info rows will have row[0] as a list, whereas the others will have a string
                if type(row[0]).__name__=='str':
                        #print "string happens"
                        non_sample_info.append(row)
                else:
                        sample_info.append(row[0])
        #print "non_sample_info before", non_sample_info
        #print "sample_info", sample_info                        

        #prune non_sample_info of 'Pop' or equivalent

        non_sample_info=filter (lambda a: a not in ['Pop','POP','pop', ['Pop'],['POP'],['pop']] , non_sample_info)
        #print "pruned", non_sample_info


        #prune non_sample_info of population_names

        #print "population_names", population_names
        non_sample_info=filter(lambda a: a[0] not in population_names, non_sample_info)
        #print "after removing population names", non_sample_info

        #combine the loci into one list
        header=[row[0] for row in non_sample_info]
        #print "header", header

        #insert 'ID' into front of header
        header.insert(0,'ID')
        #print "header", header


        #convert all -99 values in sample_info into 'NA'
        for row in sample_info:
                for i, o in enumerate(row):
                        if o==-99:
                                row[i]='NA'
                                

        #remove the commas from the first element in each row of sample_info
        for row in sample_info:
                row[0]=row[0].replace(',','')
        #print "commas removed", sample_info

        #add header to sample_info
        sample_info.insert(0, header)
        #print "with header", sample_info

        #save the output
        sample_info==remove_spaces_from_header(sample_info)
        saveCsv(sample_info, 'heterozygotes_and_homozygotes.csv', pipeline_directory)

def Incorporate_MLH_into_everythingFile(pipeline_directory):
        donor=readCsv_input(pipeline_directory+'MLH_calculated.csv')
        #print "donor",donor, "done now"
        acceptor=readCsv_input(pipeline_directory+'MLH_output.csv')
        #print "acceptor", acceptor

        acceptorTrans=transposed(acceptor)
        #print "acceptorTrans", acceptorTrans

        #print donor
        tally=0
        for i in range (1,len(donor)):
            if donor[i][0] in acceptorTrans[0]:
                tally+=1
                acceptor[acceptorTrans[0].index(donor[i][0])].append(donor[i][1])#.replace(" ",'')
                #acceptor[acceptorTrans[0].index(donor[i][0])].append(donor[i][2])
                #acceptor[acceptorTrans[0].index(donor[i][0])].append(donor[i][3])
        #print acceptor
        #print "tally", tally

        #add IR, SH, HL columns
        acceptor[0].append('MLH')
        #acceptor[0].append('SH')
        #acceptor[0].append('HL')
        #print "acceptor before", acceptor
        acceptor=fill_blanks_of_data_table(acceptor)
        acceptor=make_sure_NA_are_numeric_compliant(acceptor)
        acceptor=remove_spaces_from_header(acceptor)
        saveCsv(acceptor,"MLH_output.csv", pipeline_directory)

        print "Now just plain MLH values have also been incorporated into the fitness data file (called MLH_output.csv)"

def make_sure_NA_are_numeric_compliant(data_original):
        data=copy.deepcopy(data_original)
        for row in data:
                for index, item in enumerate(row):
                        if item ==' NA':
                                row[index]='NA'
        return data

def generate_trait_spreadsheet (directory_address_with_csv_files):
    master_data=[['ID','Reg_coeff']]
    file_list=os.listdir(directory_address_with_csv_files)
    only_csvs=[x for x in file_list if '.csv' in x]
    
    #assumes you want the second row, and that you want to strip everything that comes before the tab/several spaces
    for row in only_csvs:
        #print "row", row
        sample=[x for x in row if x.isdigit()]
        sample_str=''.join(sample)
        #print "sample_str", sample_str
        data=readCsv_input(directory_address_with_csv_files+row)
        #del data[0]
        print data[1]
        coef_candidates=data[1][0].split()
        regression_coef=coef_candidates[len(coef_candidates)-1]
        #print regression_coef
        master_data.append([sample_str,regression_coef])

    saveCsv(master_data,'final.csv',directory_address_with_csv_files)

def remove_spaces_from_all_entries(list_of_lists):
        col_list=[]
        row_list=[]
        for row in list_of_lists:
                #print "row", row
                for entry in row:
                        #print "entry", entry
                        row_list.append(entry.replace(" ",""))
                col_list.append(row_list)
                row_list=[]
                
        #print copy_list
        return col_list

def combine_het_homo_and_fitness_scores(address_of_het_homo, address_of_fitness, pipeline_directory):
 #       het_homo=readCsv_input(address_of_het_homo)
#        fitness=readCsv_input(address_of_fitness)
        incorporate_donor_into_acceptor(address_of_het_homo, address_of_fitness,pipeline_directory+'genotype_fitness_combined_for_slh_test.txt')

        #retrieve that text file
        lines=open(pipeline_directory+'genotype_fitness_combined_for_slh_test.txt').readlines()
        #print lines[0]
        new_lines=[]

        #make it into a csv file
        for line in copy.deepcopy(lines):
                #print "line before", line
                line=line.replace("'","")
                #print "line after", line
                line=line.replace("[","")
                line=line.replace("]","")
                line=line.replace("\n","")
                #print line
                new_lines.append(line)
                
        #print "before split", new_lines
        new_lines=[x.split(',') for x in new_lines]
        #print "after split", new_lines
        new_lines=remove_spaces_from_header(new_lines)
        all_spaces_removed=remove_spaces_from_all_entries(new_lines)
        saveCsv(all_spaces_removed, 'genotype_fitness_combine_for_slh_test_as_spreadsheet.csv',pipeline_directory)

def replace_internal_space(test_string):
    try:
        matched=re.search("\w+\s\w+", test_string)
        return "_".join(matched.group(0).split())
    except:
        return test_string

def run_model_comparison(csv_address, keeplist_address, pipeline_directory):
        try:
                os.mkdir(pipeline_directory+'Single_Locus_Effects_Test')
        except:
                pass

        r=R(use_numpy=True)

        print r("library(pgirmess)")
        print r("library(MASS)")
        
        #convert NAs in locus-containing columns into averages

        keeplist=readCsv_input(keeplist_address)
        #print "keeplist", keeplist
        keeplist_as_list_of_strings=[x[0].replace(' ','_') for x in keeplist]
        #print "keeplist_as_list_of_strings", keeplist_as_list_of_strings
        keeplist_as_list_of_strings=[x.replace('-','_') for x in keeplist_as_list_of_strings]
        #re.sub(r'(\w)-(\w)', lambda m: '_' (m.groups()), x)
        #print "keeplist_as_list_of_strings", keeplist_as_list_of_strings

        master_file=readCsv_input(csv_address)
        master_file=remove_spaces_from_all_entries(master_file)
        #master_file=remove_spaces_from_header(master_file)
        #print "master_file", master_file
        master_trans=transposed(copy.deepcopy(master_file))
        #print "master_trans", master_trans
        #print "locus names:\n"
        for x in master_trans:
                #print "x[0] before", x[0]
                x[0]=replace_internal_space(x[0])
                #print "x[0] after internal space replace", x[0]
                x[0]=x[0].replace('-','_')
                x[0]=x[0].replace(" ", "")
                #print "x[0] after both", x[0]
        #print "master_trans", master_trans
        just_loci=[x for x in master_trans if x[0].replace(' ','') in keeplist_as_list_of_strings]
        just_fitness=[x for x in master_trans if x[0].replace(' ','') not in keeplist_as_list_of_strings]
        just_fitness_names=[x[0] for x in just_fitness]
        #print "just_fitness", just_fitness
        #print "just_loci", just_loci
        #print "just_fitness_names", just_fitness_names
        #print "just_fitness_names", just_fitness_names
        NAs_filled=[score_NAs(row) for row in just_loci]
        #print "NAs_filled", NAs_filled
        #NAs_filled_flipped_back=transposed(NAs_filled)
        #print NAs_filled_flipped_back
        #saveCsv(NAs_filled_flipped_back, 'hets_and_homos_with_NAs_filled_as_means.csv', pipeline_directory)

        #now replace the unfilled columns of master_file_trans with the corresponding filled ones of NAs_filled
        for i, row in enumerate(NAs_filled):
                for index, x in enumerate(master_trans):
                        ####print "row[0]", row[0], "x[0]", x[0]
                        if x[0]==row[0]:
                                master_trans[index]=NAs_filled[i]
        ready_for_comparison=transposed(master_trans)
        #print "ready_for_comparison", ready_for_comparison
        saveCsv(ready_for_comparison, 'fitness_and_loci_combined.csv',pipeline_directory+'Single_Locus_Effects_Test/')


        #the R function:
        funct="""
        make_and_compare_models <- function(fitness_trait_name, data_frame_name, vector_for_multiple_regression, predictor_for_single_regression, fam){
	#print(data_frame_name)
	#str(data_frame_name)
	sink(file=paste(paste("Multiple_regression_vs_single_regression_",fitness_trait_name, sep=''), ".txt", sep=''))
	print (fitness_trait_name)
	#this is the multiple regression and will be very specific
	fit1<-glm(data=data_frame_name, formula=as.formula(paste(fitness_trait_name,"~", paste(vector_for_multiple_regression, collapse="+"))), family=fam)
	print ("summary fit 1")
	print(summary(fit1))
	fit2<- glm(data=data_frame_name, formula=as.formula(paste(fitness_trait_name,"~",predictor_for_single_regression)), family=fam)

	print("summary fit 2")
	print(summary(fit2))
	print("model comparison stats:")
	mod_test<-anova(fit2,fit1)
	print(anova(fit2,fit1))
	print ("significance:")
        print (1-pchisq( abs(mod_test$Deviance[2]),abs(mod_test$Df[2])))
	sink()
	
	
	tiff(filename=paste(paste("Multiple_regression_",fitness_trait_name, sep=''), ".tiff", sep=''))
	layout(matrix(c(1,2,3,4),2,2))
	plot(fit1)
	dev.off()
	
	tiff(filename=paste(paste("Single_regression_",fitness_trait_name, sep=''), ".tiff", sep=''))
	layout(matrix(c(1,2,3,4),2,2))
	plot(fit2)
	dev.off()	
        }
        """


#	print (mod_test$Deviance[2])
#	print (mod_test$Df[2])
#	print (pchisq(abs(mod_test$Deviance[2]),abs(mod_test$Df[2])))


        print (r(funct))

        #implementing all of these preparations
        print "You're about to run a test that determines whether a there is evidence for a single locus effect on a particular fitness-associated trait.print. The traits:", just_fitness_names

 

        #get a list of all of the fitness traits from the user that they want to do the slh significance test on
        print "Please enter the name of each trait you'd like to run this test on. Hit return between each entry. Type 'done' when you're done. (BE MINDFUL OF SPACES!)"
        trait_list=[]
        new_value='not done'
        while new_value!='done' and new_value!='Done':
                #print "new_value", new_value, "type", type(new_value)
                new_value=raw_input()
                if (new_value not in just_fitness_names and new_value not in ['done','Done']):
                        print "You have entered a trait that is not found in the trait-containing spreadsheet. Please re-enter, being mindful of spaces and typos."
                else:
                        trait_list.append(new_value) # changed from traits_and_values.append(split_new_value_remove_spaces)
        trait_list=trait_list[0:len(trait_list)-1]
        print "trait_list", trait_list

        #get the MLH index
        print "Please enter the multilocus heterozygosity estimator you'd like to use (MLH is recommended), being mindful of spaces."
        new_value='not done'
        new_value=raw_input()
        while new_value not in just_fitness_names:
                print "You have entered a trait that is not found in the trait-containing spreadsheet. Please re-enter, being mindful of spaces and typos."
                new_value=raw_input()
        print "MLH index is", new_value
        r.mlh_value=new_value

        r.pipeline_dir=pipeline_directory
        print r("setwd(paste(pipeline_dir,'Single_Locus_Effects_Test/', sep=''))")
        print r("data<-read.csv(paste(pipeline_dir,'Single_Locus_Effects_Test/fitness_and_loci_combined.csv', sep=''),header=T)")
        r.the_loci=keeplist_as_list_of_strings
        #print r("print(data)")
        for trait in trait_list:
                #print "trait", trait
                #print type(trait)
                r.current_trait=trait

                #get the relevent link function for the response variable (and predictor variables?)
                print "The trait currently being tested is: ",trait
                family=raw_input("What sort of error distribution and link function do you expect to use in your models? Your choices are: binomial, gaussian, Gamma, inverse.gaussian, poisson, quasi, quasibinomial, and quasipoisson. Note that gaussian denotes normal linear regression. Please type one of them and press enter")
                while family not in ['binomial', 'gaussian', 'Gamma', 'inverse.gaussian', 'poisson', 'quasi', 'quasibinomial', 'quasipoisson']:
                        family=raw_input("You typed a link function that is not one of the options. Please type either binomial, gaussian, Gamma, inverse.gaussian, poisson, quasi, quasibinomial, or quasipoisson. Note that gaussian denotes normal linear regression.")
                r.fam=family
                #print r("print(the_loci)")
                #print r("print(data)")
                #print r("print(current_trait)")
                #print r("print(data[current_trait])")
                print (r('make_and_compare_models(current_trait, data[!is.na(data[current_trait]),], the_loci, mlh_value, fam)'))
                                              

def score_NAs(the_list): #where the_list is a list such that the first element is a locus name and the other elements are the heterozygosity status of each sample at that locus. This module fills the NAs with the average heterozygosity of those scored
        list_copy=copy.deepcopy(the_list)
        #print the_list[0]
        del list_copy[0]
        #print list_copy

        #remove spaces and make ints
        as_ints=[]
        for x in list_copy:
                x=x.replace(' ','')
                if x in ['1', '0']:
                        x=int(x)
                        as_ints.append(x)
                else:
                        #print x, " is not an integer"
                        as_ints.append(x)
        #print "as_ints", as_ints
        total_typed=0
        total_het=0
        for row in as_ints:
                if row ==1:
                        total_typed+=1
                        total_het+=1
                if row ==0:
                        total_typed+=1
        average_for_locus=float(total_het)/float(total_typed)
        converted=[]
        for row in as_ints:
                if row in [1,0]:
                        converted.append(row)
                else:
                       converted.append(average_for_locus)
        #print "converted before", converted
        converted=insert(converted,[the_list[0]],0)
        #print "converted after", converted
        
        return converted

def all_correlations_and_pvalues(MLH_output_csv, save_directory, rejected_samples_address):
        r=R(use_numpy=True)
        r.fitness_directory=MLH_output_csv
        r.save=save_directory
        f1="""
        cor.prob<-function (X, dfr=nrow(X)-2){
	R<-cor(X, use="pairwise")
	above<-row(R) < col(R)
	r2<-R[above]^2
	Fstat<-r2*dfr/(1-r2)
	R[above]<-1-pf(Fstat,1,dfr)
	R[row(R)==col(R)]<-NA
	R
        }
        """
        f2="""
        cor.prob.spearman<-function (X, dfr=nrow(X)-2){
	R<-cor(X, method='spearman', use="pairwise")
	above<-row(R) < col(R)
	r2<-R[above]^2
	Fstat<-r2*dfr/(1-r2)
	R[above]<-1-pf(Fstat,1,dfr)
	R[row(R)==col(R)]<-NA
	R
        }
        """
        f3="""
        flattenSquareMatrix<-function(m){
	if((class(m) !='matrix')|(nrow(m) != ncol(m))) stop("Must be a square matrix.")
	if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
	ut<-upper.tri(m)
	data.frame(i=rownames(m)[row(m)[ut]],
	j=rownames(m)[col(m)[ut]],
	cor=t(m)[ut],
	p=m[ut])
        }
        """
        f4="""
        cor.test.p <- function(x){
            FUN <- function(x, y) cor.test(x, y)[[3]]
            z <- outer(
              colnames(x), 
              colnames(x), 
              Vectorize(function(i,j) FUN(x[,i], x[,j]))
            )
            dimnames(z) <- list(colnames(x), colnames(x))
            z
        }
        """
        f5="""
        cor.test.p.s <- function(x){
            FUN <- function(x, y) cor.test(x, y, method='spearman')[[3]]
            z <- outer(
              colnames(x), 
              colnames(x), 
              Vectorize(function(i,j) FUN(x[,i], x[,j]))
            )
            dimnames(z) <- list(colnames(x), colnames(x))
            z
        }
        """
 
        print(r(f1))
        print(r(f2))
        print(r(f3))
        print(r(f4))
        print r(f5)
        print r("data_for_cor<-read.csv(fitness_directory, header=T)")#this could be a problem line
        r.rejected_path=rejected_samples_address
        print r("rejected<-read.csv(rejected_path, header=F)")
        print r("data_for_cor<-subset(data_for_cor, !(ID %in% rejected[,1]))") 
        print r("just_nums<-sapply(data_for_cor,is.numeric)")
        print r("just_nums_data<-data_for_cor[,just_nums]")
        #print r("print ('just_nums_data')")
        print r("print (just_nums_data)")
        print r("print (str(just_nums_data))")
        print r("z<-cor.test.p(just_nums_data)")
        print r("w<-cor.prob(just_nums_data)")
        #print r("print(w)")
        print r("x<-z")
        print r("x[lower.tri(x)]<-w[lower.tri(w)]")
        
        print r("x<-flattenSquareMatrix(x)")
        #print r("print(x)")
        print r("x=data.frame(x,'pearson_adjusted_p'=p.adjust(x$p, 'BH'))")
        print r("names(x)[3]<-'pearson_cor'")
        print r("names(x)[4]<-'pearson_p'")


        print r("a<-cor.test.p.s(just_nums_data)")
        print r("b<-cor.prob.spearman(just_nums_data)")
        print r("c<-a")
        print r("c[lower.tri(c)]<-b[lower.tri(b)]")
        print r("y<-flattenSquareMatrix(c)")
        #print r("print(y)")
        
        #print r("y<-flattenSquareMatrix(cor.test.p.s(just_nums_data))")
        print r("y=data.frame(y,'spearman_adjusted_p'=p.adjust(y$p, 'BH'))")
        print r("names(y)[3]<-'spearman_cor'")
        print r("names(y)[4]<-'spearman_p'")
        print r("z<-x[,1:2]")
        print r("z$pearson_corr<-x[,3]")
        print r("z$spearman_corr<-y[,3]")
        print r("z$pearson_p<-x[,4]")
        print r("z$pearson_p_adjusted<-x[,5]")
        print r("z$spearman_p<-y[,4]")
        print r("z$spearman_p_adjusted<-y[,5]")
        print r("write.csv(z,file=paste(save,'all_correlations.csv',sep=''))") #this could be a problem line
                
def generate_allele_frequency_spreadsheet(genePop_file_address, pipeline_directory):
    linelist=open(genePop_file_address).readlines()
    #print type(linelist)
    #print linelist
    #print "module 1 entered"
    good_data=[]
    locus_names=[]
    parsing = False
    new=[]
    regex_query=re.compile('Locus:\s\w+')

    ##get a list of locus names##
    for line in linelist:
            if "Locus" in line:
                    #print "this line",line
                    locus_names.append(re.search(regex_query,line).group())
    #print "after regex", locus_names
    #for line in locus_names:
#            line.split('Locus:')
#    for line in locus_names:
 #           new.append(line.replace('\n',''))
#    locus_names=[]
#    for line in new:
#            locus_names.append(line.replace('Pop: 652',''))
#    new=[]
    for line in locus_names:
            new.append(line.replace('Locus:',''))
    locus_names=remove_spaces_from_list_recursive(new)
    #print locus_names
    
    ##isolate the lines between Fis and Tot and store them in good_data##
    for line in linelist:
            if line.startswith("                                           Fis"):
                parsing = True
            elif line.startswith("    Tot"):
                parsing = False
            if parsing:
                good_data.append(line.replace('\n','')) 
    #print "good_data", good_data

    ##split the data by spacing##
    for i, entry in enumerate(good_data):
            good_data[i]=entry.split()
    #print good_data

    header=['Allele', 'Sample', 'count', 'Frequency', 'W&C', 'R&H']
    loc_count=0
    for i,entry in enumerate(good_data):
            if entry==header: 
                    good_data[i]=[locus_names[loc_count]]
                    loc_count+=1
    #print good_data

    ##remove Fis and ['----------------'] lines##
    header2=['Allele_Name', 'Sample_count', 'Allele_Frequency', 'Fis_WC']
    Fis_removed=[header2]
    for entry in good_data:
            if entry not in [['----------------'],['Fis']]:
                    Fis_removed.append(entry)

    ##insert relevant locus name in front of every row##                   
    locus_tracker='Locus_name'                   
    for x in range (0, len(Fis_removed)):
            if Fis_removed[x][0] in locus_names:
                    locus_tracker=Fis_removed[x][0]
            Fis_removed[x].insert(0,locus_tracker)
            
    #print Fis_removed

    ##remove the rows with two entries in lcous_names##
    Locus_dups_removed=[]
    for entry in Fis_removed:
            if not contains_dupes(entry):
                    Locus_dups_removed.append(entry)
    #print "gets to the end of the generate allele_spreadsheet module", Locus_dups_removed                   
    saveCsv(Locus_dups_removed,'allele_freqs.csv',pipeline_directory)

def contains_dupes(a_list):
        unique_list=[]
        for item in a_list:
                if item not in unique_list:
                        unique_list.append(item)
        return (len(unique_list)<len(a_list))

def uniques(seq): 
   # order preserving
   checked = []
   for e in seq:
       if e not in checked:
           checked.append(e)
   return checked

def effective_allele_number(allele_freqs_file_address, pipeline_directory):
   #print "script gets to beginning of this module"
   original=readCsv_input(allele_freqs_file_address)
   #print data

   data=copy.deepcopy(original)
   del data[0]
   #print data

   locus_list=[x[0] for x in data]
   locus_list=uniques(locus_list)
   #print locus_list

   effective_alleles=[['Locus','Effect_number_alleles', 'Actual_number_alleles']]
   actual_number_allles=[["Actual_number_alleles"]]
   for locus in locus_list:
       sum_p=0
       allele_count=0
       for row in data:
           #print "row[0]", row[0], "locus", locus
           if row[0]==locus:
               #print "row[0]", row[0], "locus", locus
               sum_p+=(float(row[3])*float(row[3]))
               allele_count+=1
       #print "sum_p", sum_p
       #print "allele_count", allele_count
       effective_alleles.append([locus, float(1.0/sum_p), allele_count])
       #actual_number_alleles.append(allele_count)

   #data_trans=HefPipe_modules.transposed(original)
   #data_trans.append(effective_alleles)
   #back=HefPipe_modules.transposed(data_trans)
   #print "script gets here", effective_alleles
   saveCsv(effective_alleles, 'Effective_alleles_per_locus.csv', pipeline_directory)                               

def generate_obs_exp_het_homo_spreadsheet(genePop_file_address, pipeline_directory):
    linelist=open(genePop_file_address).readlines()
    #print type(linelist)
    #print linelist
    #print "module 1 entered"
    good_data=[]
    locus_names=[]
    parsing = False
    new=[]
    regex_query=re.compile('Locus:\s\w+')

    ##get a list of locus names##
    for line in linelist:
            if "Locus" in line:
                    #print "this line",line
                    locus_names.append(re.search(regex_query,line).group())
    #print "after regex", locus_names
    #for line in locus_names:
#            line.split('Locus:')
#    for line in locus_names:
 #           new.append(line.replace('\n',''))
#    locus_names=[]
#    for line in new:
#            locus_names.append(line.replace('Pop: 652',''))
#    new=[]
    for line in locus_names:
            new.append(line.replace('Locus:',''))
    locus_names=remove_spaces_from_list_recursive(new)

    locus_names=uniques(locus_names)
    
    ##isolate the four expected number of lines##
    for line in linelist:
            if line.startswith("    Expected number of homozygotes"):
                good_data.append(line.replace('\n',''))
            elif line.startswith("    Observed number of homozygotes"):
                good_data.append(line.replace('\n',''))
            elif line.startswith("    Expected number of heterozygotes"):
                good_data.append(line.replace('\n',''))
            elif line.startswith("    Observed number of heterozygotes"):
                good_data.append(line.replace('\n',''))
    #print "good_data obs het homo", good_data

    ##split the data by the colon##
    for i, entry in enumerate(good_data):
            good_data[i]=entry.split(':')
    #print "here we are",good_data

    ##remove spaces from good_data##
    good_data=remove_spaces_from_list_recursive(good_data)
    #print "good_data",good_data

    ##go through the good_data, four elements at a time, calculate H exp and Hobs, remove the four elements, and repeat for each entry in locus list##
    good_data_copy=copy.deepcopy(good_data)
    #print "perhaps locus_names has changed?", locus_names

    ##make each item in locus_names a list##
    locus_names=[[x]for x in locus_names]
    
    for locus in locus_names:
            #print "locus", locus
            H_obs=float(good_data_copy[3][1])/(float(good_data_copy[1][1])+float(good_data_copy[3][1]))
            #print "H_obs", H_obs
            locus.append(H_obs)
            H_exp=float(good_data_copy[2][1])/(float(good_data_copy[2][1])+float(good_data_copy[0][1]))
            #print "H_exp",H_exp
            locus.append(H_exp)
            del good_data_copy[0]
            del good_data_copy[0]
            del good_data_copy[0]
            del good_data_copy[0]
            #print len(good_data_copy)
    #print "new locus_names", locus_names
    locus_names.insert(0,['Locus','H_obs','H_exp'])
    #print locus_names
    saveCsv(locus_names,'H_obs_and_H_exp.csv',pipeline_directory)

def process_all_fullModel_files_in_directory(directory_address):
        dir_list=os.listdir(directory_address)
        #print "before", dir_list
        for i, f in enumerate(dir_list):
                #print f
                os.rename(directory_address+f,directory_address+f.replace(" ",""))
        dir_list=os.listdir(directory_address)
        #print type(dir_list)
        #print len(dir_list)
        fullModels=[x for x in dir_list if "fullModel.txt" in x]
        #print fullModels
        try:
                os.mkdir(directory_address+'processed_model_output/')
        except:
                pass

        for filename in fullModels:
                saveCsv(process_single_fullModel_file(directory_address+filename), title_from_address(directory_address+filename)+'_processed.csv' ,directory_address+'processed_model_output/')

def process_single_fullModel_file(file_address):
        #print "file_address", file_address
        f=open(file_address).readlines()
        the_cream=[]
        parse=False
        for line in f:
                #print "line", line
                if line.startswith('Coefficients'):
                        parse=True
                        #print "starts with Coefficients", line
                if line.startswith('--'):
                        parse=False
                        #print "starts with --", line
                elif line.startswith('(Dispersion'):
                        parse=False
                        #print "starts with (Dispersion", line
                if parse==True:
                        #print "gets added:", line
                        the_cream.append(line)
        the_cream=[x.split() for x in the_cream]
        #print "the_cream after removal of everything besides the significance details", the_cream


        ##get rid of any completely empty elements
        the_cream=[x for x in the_cream if len(x)>0]

        ##remove lines until '(Intercept)' happens
        the_cream_2=copy.deepcopy(the_cream)
        intercept_count=0
        #print "right before (Intercept) removal", the_cream_2
        for line in the_cream:
                if line[0]=='(Intercept)':
                        intercept_count+=1
                if intercept_count<1:
                        del the_cream_2[0]

        #print "right before the problem", the_cream_2
        ##'(Intercept)' should be the first element of the first line now; delete this line, too
        del the_cream_2[0]
        #print the_cream_2               
        
        the_cream_2=remove_spaces_from_list_recursive(the_cream_2)

        ##remove the significance values
        for line in the_cream_2:
                for entry in line:
                        if entry in ['*','**', '***','.', '<']:
                                line.remove(entry)
                                


        ##append the type of score (z or t, for example) to the end of each row
        score_type=the_cream[1][3]
        #print "score_type", score_type
        for line in the_cream_2:
                line.append(score_type)

        ##alert user of errors
        for line in the_cream_2:
                if len(line)!=6:
                        print "something is wrong with the glm or lm file being processed", the_cream_2
        
        ##append header row
        the_cream_2.insert(0,['Trait','Estimate','Std_error','score','p_value', 'type_of_score'])#score can be a z-value or a t-value
        return the_cream_2

def title_from_address(string_of_address):
        #print "string_of_address", string_of_address
        query=re.compile(r'/(\w*)\.txt')
        search1=re.search(query,string_of_address)
        return search1.group(1)

def generate_gephast_files(fitness_file_address, genotype_file_address, directory_address, rejected_samples_address):
        
        try:
                os.mkdir(directory_address+'GEPHAST_ready/')
        except:
                pass
        exceptions=readCsv_input(rejected_samples_address)

        ##for every column but the first (IDs) in the fitness data file
        fitness_data=readCsv_input(fitness_file_address)
        fitness_transposed=transposed(fitness_data)
        #del fitness_transposed[0] ##don't do that just yet
        #print "fitness_transposed", fitness_transposed
        for column in fitness_transposed:
                #print "column", column
                make_gephast_file([[fitness_transposed[0][x],column[x]] for x in range(0,len(fitness_transposed[0]))],genotype_file_address, directory_address+'GEPHAST_ready/', exceptions)


def make_gephast_file(list_of_trait_values, genotype_file_address, directory_address, exceptions_table): ##assumes one population at the moment
        #print"list_of_trait_values", list_of_trait_values
        #print "genotype_file_address", genotype_file_address
        #print "directory_address", directory_address
        #print "exceptions_table", exceptions_table
        genotype_file=readCsv_input(genotype_file_address)
        #print "genotype_file", genotype_file
        header=genotype_file[0]
        #print "header", header
        gephast_header=copy.deepcopy(header)
        gephast_header.remove('ID')
        gephast_header.insert(0,list_of_trait_values[0][1])
        gephast_header.insert(0,'1') #assumed one population at the moment; hence the 1s
        #print "gephast_header", gephast_header


        ##populate the rows with 1, fitness data, then genotype data (but only for samples not in the rejected lists
        final_table=[]
        final_table.append(gephast_header)
        del list_of_trait_values[0]
        for sample in list_of_trait_values:
                if [sample[0]] not in exceptions_table:
                        #print [sample[0]], "not in exceptions"
                        for row in genotype_file:
                                if row[0] == sample[0]:#sample IDs match
                                        if sample[1]!='NA': #removes samples with missing trait values
                                                final_table.append(make_gephast_list(sample,row))
        #print "after adding", final_table
        saveCsv(final_table,gephast_header[1]+'.csv',directory_address)
        return final_table
        ##replace *s with 0s

def make_gephast_list(sample,row):#a sample id, fitness trait (both in sample) and a list of genotypes for that same sample, return a list that contains ['1','trait_value','genotype1_allele1','genotype1_allele2',etc.,etc.] with ** replaced with 0
        new_table=['1',sample[1]]
        for x in range (1,len(row)):
                if row[x] in ['*','**','***','****']:
                        row[x]='0'
                new_table.append(row[x])
        #print "new_table", new_table
        return new_table

def is_row_of_blanks(row):
    #print "row pased to row_of_blanks module", row
    for entry in row:
        #print "entry", entry
        if entry !='':
            return False
    return True

def gephast_process_p_vals(directory, pipeline_directory):
        file_names=[x for x in os.listdir(directory) if '.csv' in x]
        #print "file_names", file_names

        new_table=[]
        for f in file_names:
                data=readCsv_input(directory+f)
                #print data
                new_table.append(uniques(data[0][1:len(data[0])]))
                for row in data:
                    if row[0]=='':
                        if not is_row_of_blanks(row):               
                            #print "row", row
                            new_row=[]
                            new_row=[x for x in row if x!='']
                            #print "new_row", new_row
                            new_row.insert(0,'')
                            #print "new_row", new_row
                            #new_row.insert(0,'')
                            new_table.append(new_row)
        #print "new_table", new_table
#        saveCsv(new_table,'gephast_pvals.csv',pipeline_directory)

        ##get the header of the master file
        master_table=[]
        master_header=[]
        for row in new_table:
            if row[0]!='':
                    master_header.append(row[0])
        #print "Master_header", master_header                    

        ##for every item besides the first in row in new_table where the first entry is not '', append to it the float version of the element in the next row /100
        for i, row in enumerate(new_table):
                if row[0]!='':
                        for index,entry in enumerate(row):
                                if index==0:
                                        master_table.append([entry])
                                if index!=0:
                                        master_table.append([entry,float(new_table[i+1][index])/100])
        #print "master_table", master_table

        ##now, you have a list of lists, such that each trait name is an element (i.e. list of len 1) followed by each locus gephasted for that trait with pvalue (lists of len 2), followed by the next trait name...                
        ##make it so that each list in new_master is a list of trait name followed by each locus and pvalue (i.e., separate the trait-plus-locus/pvalues into separate lists)
        new_master=[]
        count=1
        for row in master_table:
                #print "row", row
                if len(row)==1:
                        new_master.append([[row[0],'pval'+str(count)]])
                        count+=1
                else:
                        new_master[len(new_master)-1].append(row)
        #print "new_master", new_master

        #for row 
        nm_trans=transposed(new_master)
        #print "nm_trans", nm_trans
        #saveCsv_v2(nm_trans,'gephast_pre_process_text.csv',pipeline_directory)
        nm_take_out=[]
        temp=[]
        for row in nm_trans:
                #print "row", row
                temp=list_of_lists_to_list(row)
                #print "temp", temp[0]
                nm_take_out.append(temp)
                temp=[]
        #print nm_take_out
        saveCsv_v2(nm_take_out,'gephast_pvals.csv',pipeline_directory)
        #print test

        ##and now create p-adjusted columns
        r=R(use_numpy=True)
        r.pipel_directory=pipeline_directory

        print r("data<-read.csv(paste(pipel_directory,'gephast_pvals.csv',sep=''), header=T)")
        print r("sub1<-data[c(F,T)]")
        print r("sub2<-data[c(T,F)]")
        longCMD="""

        for (i in 1:ncol(sub1)){
        	sub1[,i]=p.adjust(sub1[,i],'BH')
        }
        """
        print r(longCMD)
        print r("data[c(F,T)]<-sub1")
        print r("write.csv(data,file=paste(pipel_directory,'gephast_p_adjusted.csv',sep=''))")

def is_list_of_strings(the_list):
        for item in the_list:
                if not isinstance(item, basestring):
                        return False
                else:
                        return True

def list_of_lists_to_list(list_of_lists):
        #print "before", list_of_lists
        final=[]
        try:
                for sublist in list_of_lists:
                        try:
                                for item in sublist:
                                        final.append(item)
                        except:
                                #print "exception for ",sublist
                                final.append(sublist)
                                #has to be twice because in this case both the locus and the p-value are missing
                                final.append(sublist)
        except:
                #print "super exception for ",list_of_lists
                final.append(list_of_lists)
        #print "after", final
        return final

def spearman_corr_chart(pipeline_directory, address_of_acceptor_file, rejected_samples_address):
        r=R(use_numpy=True)
        r.pipeline_dir=pipeline_directory
        r.acceptr_address=address_of_acceptor_file
        print r("library('PerformanceAnalytics')")
        #print r("tmpstr <- deparse(chart.Correlation)")
        cmd1="""
	my.spearman.chart.correlation<-function (R, histogram = TRUE, ...){
		x = checkData(R, method = "matrix")                                                   
		panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs", 	cex.cor.scale=1, cex.cor, ...) {
			usr <- par("usr")                                                                 
	    	on.exit(par(usr))                                                                   
	    	par(usr = c(0, 1, 0, 1))                                                            
	    	r <- (cor(x, y,method='spearman', use = use))                                                      
		    txt <- format(c(r, 0.123456789), digits = digits)[1]                                
		    txt <- paste(prefix, txt, sep = "")                                              
		    if (missing(cex.cor)){
		    	cex <- 0.8/strwidth(txt)
		    }                                                                                                                                  
		        test <- cor.test(x, y, method='spearman')                                                             
		        Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))                                               
		        text(0.5, 0.5, txt, cex = cex * r^cex.cor.scale)                                                  
		        text(0.8, 0.8, Signif, cex = cex, col = 2)                                           
		}                                                                                       
		f <- function(t) {
			dnorm(t, mean = mean(x), sd = sd.xts(x))                                            
		}                                                                                       
		hist.panel = function(x, ...) {
			par(new = TRUE)                                                                     
		    hist(x, col = "light gray", probability = TRUE, axes = FALSE, main = "", breaks = "FD")                                                   
		    lines(density(x, na.rm = TRUE), col = "red", lwd = 1)                             
		    rug(x)                                                                              
		}                                                                                       
		if (histogram){
			pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = hist.panel, ...)
		}                                                   
		    else {
		    	pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, ...)
		    }                                                                                
	} 
        
        """
        print r(cmd1)
#        print r("abscor<-grep('^ *r <- abs',tmpstr)")
        #print r("(print(abscor))")
        #print r(cmd2)
#        print r("tmpstr[abscor]<-paste('r<-cor(x,y,method='spearman',use=use)')")
        #print r("panelcorline <- grep('^ *panel.cor',tmpstr)")
        #print r("tmpstr[panelcorline] <- paste(tmpstr[panelcorline],'cex.cor.scale=1,')")
        #print r("rscaleline <- grep('^ *text\\(0.5',tmpstr)")
        #print r("tmpstr[rscaleline] <- gsub('cex \\* r','cex*r^cex.cor.scale',tmpstr[rscaleline])")
        #print r("spearman.chart.Correlation <- eval(parse(text=tmpstr))")
        print r("data_for_cor<-read.csv(paste(acceptr_address), header=T)")
        r.rejected_path=rejected_samples_address
        print r("rejected<-read.csv(rejected_path, header=F)")
        print r("data_for_cor<-subset(data_for_cor, !(ID %in% rejected[,1]))") 
        print r("data_for_cor<-data_for_cor[,2:ncol(data_for_cor)]")#gets rid of ID column
        print r("just_nums<-sapply(data_for_cor,is.numeric)")
        #print r("print(just_nums)")
        print r("just_nums_data<-data_for_cor[,just_nums]")
        print r("library(PerformanceAnalytics)")
        print r("pdf(paste(pipeline_dir,'Correlations/corr_chart_spearman.pdf',sep=''))")
        print r("my.spearman.chart.correlation(just_nums_data, histogram=TRUE, pch='+', cex.cor.scale=0, cex.labels=0.2)")
        print r("dev.off()")



def pearson_corr_chart(pipeline_directory, address_of_acceptor_file, rejected_samples_address):
        r=R(use_numpy=True)
        r.pipeline_dir=pipeline_directory
        r.acceptr_address=address_of_acceptor_file
        print r("library('PerformanceAnalytics')")
        cmd1="""
	my.pearson.chart.correlation<-function (R, histogram = TRUE, ...)
	{
		x = checkData(R, method = "matrix")                                                   
		panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs", 	cex.cor.scale=1, cex.cor, ...) {
			usr <- par("usr")                                                                 
	    	on.exit(par(usr))                                                                   
	    	par(usr = c(0, 1, 0, 1))                                                            
	    	r <- (cor(x, y, use = use))                                                      
		    txt <- format(c(r, 0.123456789), digits = digits)[1]                                
		    txt <- paste(prefix, txt, sep = "")                                              
		    if (missing(cex.cor)){
		    	cex <- 0.8/strwidth(txt)
		    }                                                                                                                                  
		        test <- cor.test(x, y)                                                             
		        Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))                                               
		        text(0.5, 0.5, txt, cex = cex * r^cex.cor.scale)                                                  
		        text(0.8, 0.8, Signif, cex = cex, col = 2)                                           
		}                                                                                       
		f <- function(t) {
			dnorm(t, mean = mean(x), sd = sd.xts(x))                                            
		}                                                                                       
		hist.panel = function(x, ...) {
			par(new = TRUE)                                                                     
		    hist(x, col = "light gray", probability = TRUE, axes = FALSE, main = "", breaks = "FD")                                                   
		    lines(density(x, na.rm = TRUE), col = "red", lwd = 1)                             
		    rug(x)                                                                              
		}                                                                                       
		if (histogram){
			pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = hist.panel, ...)
		}                                                   
		    else {
		    	pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, ...)
		    }                                                                                
	} 
        
        """
        print r(cmd1)
        print r("data_for_cor<-read.csv(paste(acceptr_address), header=T)")
        r.rejected_path=rejected_samples_address
        print r("rejected<-read.csv(rejected_path, header=F)")
        print r("data_for_cor<-subset(data_for_cor, !(ID %in% rejected[,1]))") 
        print r("data_for_cor<-data_for_cor[,2:ncol(data_for_cor)]")#gets rid of ID column
        print r("just_nums<-sapply(data_for_cor,is.numeric)")
        #print r("names(just_nums)")
        print r("just_nums_data<-data_for_cor[,just_nums]")
        print r("library(PerformanceAnalytics)")
        print r("pdf(paste(pipeline_dir,'Correlations/corr_chart_pearson.pdf',sep=''))")
        print r("my.pearson.chart.correlation(just_nums_data, histogram=TRUE, pch='+', cex.cor.scale=0, cex.labels=0.2)")
        print r("dev.off()")


def homemade_correlation_matrix_pearson_adjusted(pipeline_directory):
        print "pearson adjusted entered"
        data=readCsv_input(pipeline_directory+'Correlations/all_correlations.csv')
        trans=transposed(data)
        culled=[x for x in trans if x[0] in ['i', 'j', 'pearson_p_adjusted', 'pearson_corr']]
        traits=uniques([x for x in culled[0] if x not in ['i', 'ID']])
        #doesn't capture the one trait that's only in the j column, in our case usually 'MLH'
        more_traits=uniques([x for x in culled[1] if x not in ['j','ID']])
        #print "more traits", more_traits
        #combine uniques of traits and more_traits
        for i in more_traits:
                if i not in traits:
                        traits.append(i)
                
        #print "new traits", traits
        

        ##create the table
        a=[[""]*i + [x]+[""]*(len(traits)-i-1) for i,x in enumerate(traits)]
        
        ##put the correlation coefficient in the upper diagonal and the significance in the lower diagonal
        cols=copy.deepcopy(culled)
        rows=copy.deepcopy(transposed(cols))
        row_dest=-1
        col_dest=-1
        for row in rows[1:]:
                for i,r in enumerate(a):
                        #print "row[0]", row[0]
                        #print "row[1]", row[1]
                        #print "r", r
                        if row[0] in r:
                                #print "this happens 1"
                                row_dest=i
                                #print "row_dest"
                        if row[1] in r:
                                #print "this happens 2"
                                col_dest=i
                        if row_dest>=0 and col_dest >=0:
                                #print "this happens 3"
                                a[row_dest][col_dest]=row[2]
                                a[col_dest][row_dest]=row[3]
                row_dest=-1
                col_dest=-1
        #print "a after", a
        saveCsv(a, 'correlation_chart_pearson_adjusted.csv', pipeline_directory+'Correlations/')                                                


def homemade_correlation_matrix_spearman_adjusted(pipeline_directory):
        print "spearman adjusted entered"
        data=readCsv_input(pipeline_directory+'Correlations/all_correlations.csv')
        trans=transposed(data)
        culled=[x for x in trans if x[0] in ['i', 'j', 'spearman_p_adjusted', 'spearman_corr']]
        traits=uniques([x for x in culled[0] if x not in ['i', 'ID']])
        #doesn't capture the one trait that's only in the j column, in our case usually 'MLH'
        more_traits=uniques([x for x in culled[1] if x not in ['j','ID']])
        #print "more traits", more_traits
        #combine uniques of traits and more_traits
        for i in more_traits:
                if i not in traits:
                        traits.append(i)
                
        #print "new traits", traits


        ##create the table
        a=[[""]*i + [x]+[""]*(len(traits)-i-1) for i,x in enumerate(traits)]
        
        ##put the correlation coefficient in the upper diagonal and the significance in the lower diagonal
        cols=copy.deepcopy(culled)
        rows=copy.deepcopy(transposed(cols))
        row_dest=-1
        col_dest=-1
        for row in rows[1:]:
                for i,r in enumerate(a):
                        #print "row[0]", row[0]
                        #print "row[1]", row[1]
                        #print "r", r
                        if row[0] in r:
                                #print "this happens 1"
                                row_dest=i
                                #print "row_dest"
                        if row[1] in r:
                                #print "this happens 2"
                                col_dest=i
                        if row_dest>=0 and col_dest >=0:
                                #print "this happens 3"
                                a[row_dest][col_dest]=row[2]
                                a[col_dest][row_dest]=row[3]
                row_dest=-1
                col_dest=-1
        #print "a after", a
        saveCsv(a, 'correlation_chart_spearman_adjusted.csv', pipeline_directory+'Correlations/')

def homemade_correlation_matrix_pearson(pipeline_directory):
        print "pearson entered"
        data=readCsv_input(pipeline_directory+'Correlations/all_correlations.csv')
        trans=transposed(data)
        culled=[x for x in trans if x[0] in ['i', 'j', 'pearson_p', 'pearson_corr']]
        traits=uniques([x for x in culled[0] if x not in ['i', 'ID']])
        #doesn't capture the one trait that's only in the j column, in our case usually 'MLH'
        more_traits=uniques([x for x in culled[1] if x not in ['j','ID']])
        #print "more traits", more_traits
        #combine uniques of traits and more_traits
        for i in more_traits:
                if i not in traits:
                        traits.append(i)
                
        #print "new traits", traits


        ##create the table
        a=[[""]*i + [x]+[""]*(len(traits)-i-1) for i,x in enumerate(traits)]
        
        ##put the correlation coefficient in the upper diagonal and the significance in the lower diagonal
        cols=copy.deepcopy(culled)
        rows=copy.deepcopy(transposed(cols))
        row_dest=-1
        col_dest=-1
        for row in rows[1:]:
                for i,r in enumerate(a):
                        #print "row[0]", row[0]
                        #print "row[1]", row[1]
                        #print "r", r
                        if row[0] in r:
                                #print "this happens 1"
                                row_dest=i
                                #print "row_dest"
                        if row[1] in r:
                                #print "this happens 2"
                                col_dest=i
                        if row_dest>=0 and col_dest >=0:
                                #print "this happens 3"
                                a[row_dest][col_dest]=row[2]
                                a[col_dest][row_dest]=row[3]
                row_dest=-1
                col_dest=-1
        #print "a after", a
        saveCsv(a, 'correlation_chart_pearson.csv', pipeline_directory+'Correlations/')                                                


def homemade_correlation_matrix_spearman(pipeline_directory):
        print "spearman entered"
        data=readCsv_input(pipeline_directory+'Correlations/all_correlations.csv')
        trans=transposed(data)
        culled=[x for x in trans if x[0] in ['i', 'j', 'spearman_p', 'spearman_corr']]
        traits=uniques([x for x in culled[0] if x not in ['i', 'ID']])
        #doesn't capture the one trait that's only in the j column, in our case usually 'MLH'
        more_traits=uniques([x for x in culled[1] if x not in ['j','ID']])
        #print "more traits", more_traits
        #combine uniques of traits and more_traits
        for i in more_traits:
                if i not in traits:
                        traits.append(i)
                
        #print "new traits", traits


        ##create the table
        a=[[""]*i + [x]+[""]*(len(traits)-i-1) for i,x in enumerate(traits)]
        
        ##put the correlation coefficient in the upper diagonal and the significance in the lower diagonal
        cols=copy.deepcopy(culled)
        rows=copy.deepcopy(transposed(cols))
        row_dest=-1
        col_dest=-1
        for row in rows[1:]:
                for i,r in enumerate(a):
                        #print "row[0]", row[0]
                        #print "row[1]", row[1]
                        #print "r", r
                        if row[0] in r:
                                #print "this happens 1"
                                row_dest=i
                                #print "row_dest"
                        if row[1] in r:
                                #print "this happens 2"
                                col_dest=i
                        if row_dest>=0 and col_dest >=0:
                                #print "this happens 3"
                                a[row_dest][col_dest]=row[2]
                                a[col_dest][row_dest]=row[3]
                row_dest=-1
                col_dest=-1
        #print "a after", a
        saveCsv(a, 'correlation_chart_spearman.csv', pipeline_directory+'Correlations/')

def num_samples_per_locus (pipeline_directory):
        data=readCsv_input(pipeline_directory+'heterozygotes_and_homozygotes.csv')
        data2=copy.deepcopy(transposed(data))
        del data2[0]
        #print "data2", data2
        final=[]
        for row in data2:
                count=0
                for x in row[1:len(row)]:
                        if x in ['1','0']:
                                count+=1
                final.append([row[0],count])
        saveCsv(final, 'number_of_samples_per_locus.csv',pipeline_directory)
        #print "final", final
