from pyper import *
import HefPipe_modules
import os
print "Make sure that you have no loci in your allele reports with only one column (more than 2 is ok, but the script is not yet equipped to deal with fewer than 2). If you don't do this, you will hopefully find out downstream due to indexing errors. If you're lucky. If you're not, beware!"

#process all of the input files
answer=raw_input("You have hopefully entered the paths to all of your files into a csv file called 'addresses.csv'. If you haven't, please do this now.Are you ready to move on? y/n")
if answer in ['y', 'Y', 'yes', 'Yes']:
        address_of_addresses=raw_input("What's the address of the addresses file? Example: /Users/mf/Desktop/addresses.csv")
        data=HefPipe_modules.readCsv_input(address_of_addresses)
        #print "addresses", data
        data_trans=HefPipe_modules.transposed(data)
        #print data_trans
        allele_reports_address=data_trans[1][0]
        keeplist_address=data_trans[1][1]
        monolist_address=data_trans[1][2]
#        desktop=data_trans[1][3]
        pipeline_directory=data_trans[1][3]
        rejected_samples_address=data_trans[1][4]
#        set_demarcation_file_address=data_trans[1][6]
#        address_of_sensitivity_GenePop_directory=data_trans[1][7]
        address_of_acceptor_file=data_trans[1][5]
        #print "got to here"
        #print (allele_reports_address,keeplist_address,monolist_address,desktop,pipeline_directory,rejected_samples_address,set_demarcation_file_address,address_of_msat_rhh,address_of_sensitivity_GenePop_directory,address_of_acceptor_file)
        #print "and got to here"
else:
        print "Please do this and run again"
        #now run the MAMs on the entire sample set
        done_status='y'
        while done_status in ['yes', 'y', 'Yes']:
                HefPipe_modules.perform_MAM(HefPipe_modules.readCsv_input(pipeline_directory+'MLH_output.csv'), pipeline_directory, rejected_samples_address)
                done_status=raw_input("Do you want to run a different MAM model? y/n")




run_Regressions=raw_input("Do you want to run regression tests?")
if run_Regressions in ['y','Y','yes','Yes']:
        #now run the regressions on the entire sample set
        done_status='y'
        while done_status in ['yes', 'y', 'Yes']:
                HefPipe_modules.perform_MAM(HefPipe_modules.readCsv_input(pipeline_directory+'MLH_output.csv'), pipeline_directory, rejected_samples_address)
                done_status=raw_input("Do you want to run a different MAM model? y/n")




        #then, run the sample-splitting bifurcations
        answer=raw_input("This pipeline is also equipped to run HFCs on only certain individuals bearing certain trait values (e.g., stressed individuals).  Would you like to run such an analysis (y/n)?")
        while answer==('y' or 'yes'):
            #create a subset of the data the subset data object is a list, where the first element is the data table, and the second element is the address of the data table
            subset_data=HefPipe_modules.get_bifurcation(HefPipe_modules.readCsv_input(pipeline_directory+'MLH_output.csv'), pipeline_directory, rejected_samples_address)

            #run the HFC models
            done_status='y'
            while done_status in ['yes', 'y', 'Yes']:
                    HefPipe_modules.perform_MAM(subset_data[0], subset_data[1], rejected_samples_address)
                    done_status=raw_input("Do you want to run a different MAM model? y/n")

            answer=raw_input("Should I continue with more subset HFCs?")
            
            #stressor=raw_input("What's the name of the column that designates stressed individuals?")
            #is there any way to pass this to the r script (is the bifuraction an r script? I can't even remember)
        print "Moving on from the stress subset analyses..."

        answer=raw_input("Do you want to run regressions on csv files of your own making (say, by merging other subsetted data?)")
        while answer in ['yes','Yes','y']:
                address_of_csv=raw_input("What's the address of the csv file you want to use?")
                directory_name=raw_input("What's the name of the directory you want to create to store these regressions in? Example: All_virus")
                try:
                        os.mkdir(pipeline_directory+directory_name+'/')
                except:
                        pass
                done_status='yes'
                while done_status in ['yes', 'y', 'Y', 'Yes']:
                        HefPipe_modules.perform_MAM(HefPipe_modules.readCsv_input(address_of_csv),pipeline_directory+directory_name+'/', rejected_samples_address)
                        done_status=raw_input("Do you want to run a different MAM model? y/n")
                answer=raw_input("Should I continue with more regressions on files of your own making?")

##run fullModel conversion to spreadsheets
answer=raw_input("Do you want to convert your fullModel.txt files from the regressions you ran into csv files?")
if answer in ['yes','y','Y','Yes']:
        HefPipe_modules.process_all_fullModel_files_in_directory(pipeline_directory+'Regressions/')
data_sub=raw_input("Do you want to convert your fullModel.txt files in the Data Subsets directory into csv files?")
if data_sub in ['yes','y','Y','Yes']:
        HefPipe_modules.process_all_fullModel_files_in_directory(raw_input("What's the path to the directory of fullModel.txt files you want to convert"))
user_made=raw_input("Do you want to convert your fullModel.txt files generated from the user-entered regression spreadsheet into csv files?")
while user_made in ['yes','y','Y','Yes']:
        HefPipe_modules.process_all_fullModel_files_in_directory(raw_input("What's the path to the directory of fullModel.txt files you want to convert"))
        user_made=raw_input("Continue with user-entered regression spreadsheets (if you had more than one)?")                        

#then, do HFC analyses on all panel-size bifurcations
print "Congratulations! You are done! You may want to re-run the analysis including or excluded different loci or individuals based on this first run, though!"
