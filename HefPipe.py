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
#clean up the allele reports to make ready for the combine module below
HefPipe_modules.allele_report_pipeline_directory_version(allele_reports_address, keeplist_address, monolist_address, pipeline_directory, rejected_samples_address)

answer=raw_input("Are you ready to move on? In other words, I'm giving you a chance to fix missing values, etc., if you need it. y/n")
if answer in ['y','Y','yes', 'Yes']:
        
        #combine all of the edited allele reports into one massive csv file and also convert this to GenePop format
        sample_names_list_e1_keys_as_numbers_e2=HefPipe_modules.combine_allele_report_pipeline_dict(pipeline_directory, keeplist_address, rejected_samples_address)
        #print "sample names pass back from combine_allele_report_pipeline to HFC_pipeline, to be passed to sensitivity_analysis", sample_names_list_e1_keys_as_numbers_e2
        print "Files have been converted to GenePop format, which can be formatted to STRUCTURE format using the program CONVERT (URL here)"

        #calculate locus statistics, such as effective number of alleles, actual number of alleles, He, and Hobs, and output them into spreadsheets
        answer=raw_input("Do you want to calculate the effective number of alleles from GenePop data that you have run through option 5 (Basic information, Fis, and gene diversities)? You'll have to have generated a .txt  file from the HTML output that you got from submitting the genePop file to GenePop Web Option 5.") 
        if answer in ['y','Y','yes', 'Yes']:
                opt5=raw_input("What's the address of the genePop option 5 file?")
                HefPipe_modules.generate_allele_frequency_spreadsheet(opt5,pipeline_directory)
                HefPipe_modules.effective_allele_number(pipeline_directory+'allele_freqs.csv', pipeline_directory)
                HefPipe_modules.generate_obs_exp_het_homo_spreadsheet(opt5,pipeline_directory)

        #convert GenePop input to rmes input
        print "converting to g2"
        HefPipe_modules.convert_to_g2 (pipeline_directory, keeplist_address)

        #convert final_output.csv to a directory of gephast-formatted files
        print "converting to GEPhAST"
        HefPipe_modules.generate_gephast_files(address_of_acceptor_file, pipeline_directory+'final_output.csv', pipeline_directory, rejected_samples_address)
        answer=raw_input("Do you want to create a list of gephast p-values for the traits for which you have run the gephast macro on?")
        if answer in ['y','Y','yes', 'Yes']:
                pause_hack=raw_input("Please do NOT proceed until you have run the gephast macro on the relevent spreadsheets and saved the output spreadsheets (these can be the same spreadsheet with the two new rows at  the bottom added) into a directory of your own making where they can be processed.")
                HefPipe_modules.gephast_process_p_vals(raw_input("What's the address of that directory? Note: all csv files in that directory will be treated as if they are the results of gephast simulations. Make sure the address ends in a backslash!"),pipeline_directory)
                ##new modules here
                
        #convert GenePop input to rhh input
        print "converting to rhh"
        HefPipe_modules.genePop_to_rhh(pipeline_directory, rejected_samples_address)

        #run rhh
        answer=raw_input("Are you ready for run the rhh R scrip y/n?")
        if answer in ['y','Y','yes', 'Yes']:
            from pyper import *
            r=R(use_numpy=True)
            r.pipeline_dir=pipeline_directory
            print r("setwd(pipeline_dir)")
            print r("require('Rhh')")
            print r("h <- mlh('rhh_ready.txt', 'Rhh_test_output.txt', '0', 4)")
            CMDS = ["sink(file='mean_and_corr_probs.txt')", "r <- h_cor('rhh_ready.txt', '0', 250, 'hl')", "mean(r)", "the_mean<-mean(r)", "print(the_mean)", "quantile(r, probs=c(0.025, 0.975))", "the_correlation<-quantile(r,probs=c(0.025, 0.975))", "print(the_correlation)", "sink()"]
            print r(CMDS)
            CMDS2= ["pdf('hhc_plot.pdf')", "hist(r)", "dev.off()"]
            print r(CMDS2)

            #incorporate the heterozygosity measures into the the master fitness-trait-containing spreadsheet
            HefPipe_modules.Incorporate_HL_from_rhh(pipeline_directory, address_of_acceptor_file)

            #and incorporate MLH (the fourth heterozygosity measure-the one not included in rhh) here
            HefPipe_modules.Incorporate_MLH_into_everythingFile(pipeline_directory)
            
        print "Moving on from MLH calculations and HHC test simulations..."


        ##generate number of samples per locus
        HefPipe_modules.num_samples_per_locus(pipeline_directory)
        
        ##correlate the fitness data

        #r=R(use_numpy=True)
        #r.address_of_acceptor_file_r=address_of_acceptor_file
        #print r("data_for_cor<-read.csv(address_of_acceptor_file_r)")
        HefPipe_modules.Correlation_heat_map(pipeline_directory, rejected_samples_address)
        HefPipe_modules.all_correlations_and_pvalues(pipeline_directory+'MLH_output.csv',pipeline_directory+'Correlations/', rejected_samples_address)
        HefPipe_modules.spearman_corr_chart(pipeline_directory, address_of_acceptor_file, rejected_samples_address)
        HefPipe_modules.pearson_corr_chart(pipeline_directory, address_of_acceptor_file, rejected_samples_address)
        HefPipe_modules.homemade_correlation_matrix_pearson_adjusted(pipeline_directory)
        HefPipe_modules.homemade_correlation_matrix_spearman_adjusted(pipeline_directory)
        HefPipe_modules.homemade_correlation_matrix_pearson(pipeline_directory)
        HefPipe_modules.homemade_correlation_matrix_spearman(pipeline_directory)
        
        #HefPipe_modules.spearman_corr_chart_padjust(pipeline_directory, address_of_acceptor_file)
        #HefPipe_modules.pearson_corr_chart_padjust(pipeline_directory, address_of_acceptor_file)


        run_slh_test=raw_input("Do you want to test for single-locus effects?")
        if run_slh_test in ['y','Y','yes', 'Yes']:
                #run the single-locus association test (f-ratio test of multiple regression model and MLH-as-single-predictor model)
                HefPipe_modules.combine_het_homo_and_fitness_scores(pipeline_directory+'heterozygotes_and_homozygotes.csv', pipeline_directory+'MLH_output.csv', pipeline_directory)
                HefPipe_modules.run_model_comparison(pipeline_directory+'genotype_fitness_combine_for_slh_test_as_spreadsheet.csv', keeplist_address, pipeline_directory)

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
                        HefPipe_modules.process_all_fullModel_files_in_directory(pipeline_directory+'/Data_subsets/')
                user_made=raw_input("Do you want to convert your fullModel.txt files generated from the user-entered regression spreadsheet into csv files?")
                while user_made in ['yes','y','Y','Yes']:
                        HefPipe_modules.process_all_fullModel_files_in_directory(raw_input("What's the path to the directory of fullModel.txt files you want to convert"))
                        user_made=raw_input("Continue with user-entered regression spreadsheets (if you had more than one)?")                        

        #then, do HFC analyses on all panel-size bifurcations
        print "Congratulations! You are done! You may want to re-run the analysis including or excluded different loci or individuals based on this first run, though!"

        ##exclude this from the actual HefPipe pub-ready script:
        import supplement
