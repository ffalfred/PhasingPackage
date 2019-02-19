import os
import math
import numpy as np
import pandas as pd
import operator
pd.options.mode.chained_assignment = None

class phasing_class_tool:

    def __init__(self,SNPSfile,min_SNPS,output_dir):
        
        """Instants of the files. Include the input file and the number of consecutive 
        SNPs required to define an haploblock"""
        
        self.inputfile = pd.read_csv(SNPSfile,sep='\t')   # Input file in .txt format
        self.min_snps = min_SNPS    # Minimum of SNPs
        self.output_dir = str(output_dir)

    def limit_snps(self,array_,pos_array):
        
        """ Function that controls the length of the haploblock. Regarding the Troubleshooting 
        Table, Step 44: Very closely positioned crossovers that give rise tovery short haplotype
        blocks. It deletes the haploblcoks shorter than the minimum amount of SNPs, and joins the
        flanking haploblocks. It only works with egg and PB2 (without heterozygosity)."""
        
        pos_pd=pd.DataFrame(pos_array).set_index('Position')    # Dataframe with position and phase of the SNPs
        array_=np.append(array_,np.zeros((1,np.shape(array_)[1])),axis=0)   # Add column to the phases array
        
        for i in  range(np.shape(array_)[1]):
            num_snps=pd.notna(pos_pd.loc[array_[1,i]:array_[2,i]][pos_pd.columns[0]]).sum() # Count amount of SNPs with defined phase
            array_[3,i]=num_snps
        
        haplo_count=pd.DataFrame(array_.T,columns=['Haplotype','Start','End','NumSNPS'])    # Define dataframe with phase, start and end of haploblock and number of defined SNPs contained
        del_haplo_count=haplo_count[(haplo_count['Haplotype']!=0.5)&(haplo_count['NumSNPS']>=self.min_snps)].reset_index()  # Delete haploblocks smaller than min_snps

        deleting_rows=[]
        # In order to join haploblocks after deleting small haploblocks, iterate over rows, comparing the row with the previous one.
        for index, row in del_haplo_count.iterrows():
            if index!=0:
                if del_haplo_count.loc[index]['Haplotype']==del_haplo_count.loc[index-1]['Haplotype']:  #If haploblocks have same phase
                    # Form an haploblock from both haploblocks
                    del_haplo_count.at[index,'Start']=del_haplo_count.loc[index-1]['Start']
                    del_haplo_count.at[index,'NumSNPS']=del_haplo_count.loc[index-1]['NumSNPS']+del_haplo_count.loc[index]['NumSNPS']
                    deleting_rows.append(index-1)
        
        del_haplo_count.drop(deleting_rows, inplace=True)   # Delete haploblocks used
        del_haplo_count=del_haplo_count.reset_index()   #Reset index
        del_haplo_count.drop('level_0',axis=1,inplace=True)
        del_haplo_count.drop('index',axis=1,inplace=True)
        return(del_haplo_count)     #Return dataframe with phase, start and end of haploblocks


    def treat_PB1(self,dict_cells,mark_pb1,i):
        """ As it might contain heterogyzous SNPs, delimiting PB1 haploblocks using the rule: Each recombination should 
        always be present in two out of three of the products of meiosis at the same position."""
        ## Select all phase changing positions in PB2 and EGG
        regions=[]
        pos_pd=mark_pb1[['Position',str(i)+'.PB1.Phase']].set_index('Position')
        for d in dict_cells:
            regions.extend(dict_cells[d]['Start'].values)
            regions.extend(dict_cells[d]['End'].values)
        regions=list(set(regions))
        regions.sort()
        # Divide the PB1 chromosome in all the phase changing positions, and analyze in which phase are
        new_coord=np.zeros((len(regions)-1,2))
        list_=[]
        for n in range(len(regions)):
            if n+1<len(regions):
                if regions[n]!=regions[n+1]:    
                    num_snps=pos_pd.loc[regions[n]:regions[n+1],str(i)+'.PB1.Phase'].count()    #Count SNPs in the block
                    # Get percentage of SNPs are not in phase
                    try:
                        value_0=pos_pd.loc[regions[n]:regions[n+1],str(i)+'.PB1.Phase'].value_counts(normalize=True).loc[0]
                    except KeyError:
                        value_0=0
                    # Get percentage of SNPs are heterozygous
                    try:
                        value_05=pos_pd.loc[regions[n]:regions[n+1],str(i)+'.PB1.Phase'].value_counts(normalize=True).loc[0.5]
                    except KeyError:
                        value_05=0
                    # Get percentage of SNPs are in phase
                    try:
                        value_1=pos_pd.loc[regions[n]:regions[n+1],str(i)+'.PB1.Phase'].value_counts(normalize=True).loc[1]
                    except KeyError:
                        value_1=0
                    # Define haploblock. If the percentage of 0 or 1 are higher than 66%, homozygous. If not, heterozygous.
                    if value_0>0.66:
                        haplo=0
                    elif value_1>0.66:
                        haplo=1
                    else:
                        haplo=0.5
                    list_.append([regions[n],regions[n+1],value_0,value_05,value_1,num_snps,haplo])
        first_haplo=pd.DataFrame(list_,columns=['Start','End','0','0.5','1','SNPS','Haplotype'])
        ## Delete haploblocks smaller than min_snps
        last_list=[]
        done_rows=[]
        for index,row in first_haplo.iterrows():
            if (row['SNPS']<self.min_snps) & (index not in done_rows):  #If haploblocks is smaller than min_snps
                if index==0 or row['Haplotype']==first_haplo['Haplotype'][index+1]:     #If it is the first haploblock, or the haploblock has the same phase as the next one
                    # Calculate the percentages of the new haploblock done with the haploblock and the next one
                    value_0=((row['0']*row['SNPS'])+(first_haplo['0'][index+1]*first_haplo['SNPS'][index+1]))/(row['SNPS']+first_haplo['SNPS'][index+1])
                    value_05=((row['0.5']*row['SNPS'])+(first_haplo['0.5'][index+1]*first_haplo['SNPS'][index+1]))/(row['SNPS']+first_haplo['SNPS'][index+1])
                    value_1=((row['1']*row['SNPS'])+(first_haplo['1'][index+1]*first_haplo['SNPS'][index+1]))/(row['SNPS']+first_haplo['SNPS'][index+1])
                    # Calculate the haplotype with the same 66% rule (two thirds)
                    if value_0>0.66:
                        haplo=0
                    elif value_1>0.66:
                        haplo=1
                    else:
                        haplo=0.5
                    last_list.append([row['Start'],first_haplo['End'][index+1],value_0,value_05,value_1,row['SNPS']+first_haplo['SNPS'][index+1],haplo])
                    done_rows.append(index+1)
                    # If it is the last haploblock, or the haploblock phase is the same as the previous, or the flanking haploblocks phases are equal
                elif index==len(first_haplo.index) or row['Haplotype']==first_haplo['Haplotype'][index-1] or first_haplo['Haplotype'][index+1]==first_haplo['Haplotype'][index-1]:
                    # Notice that, if the flanking haploblocks are equal and the haploblock smaller, the haploblock will be ignored.
                    # If they are not equal, it will be deleted.
                    # Calculate the percentages of the new haploblock done with the haploblock and the previous one
                    value_0=((row['0']*row['SNPS'])+(first_haplo['0'][index-1]*first_haplo['SNPS'][index-1]))/(row['SNPS']+first_haplo['SNPS'][index-1])
                    value_05=((row['0.5']*row['SNPS'])+(first_haplo['0.5'][index+1]*first_haplo['SNPS'][index-1]))/(row['SNPS']+first_haplo['SNPS'][index-1])
                    value_1=((row['1']*row['SNPS'])+(first_haplo['1'][index-1]*first_haplo['SNPS'][index-1]))/(row['SNPS']+first_haplo['SNPS'][index-1])
                    # Calculate the haplotype with the same 66% rule (two thirds)
                    if value_0>0.66:
                        haplo=0
                    elif value_1>0.66:
                        haplo=1
                    else:
                        haplo=0.5
                    last_list.append([first_haplo['Start'][index-1],row['End'],value_0,value_05,value_1,row['SNPS']+first_haplo['SNPS'][index-1],haplo])
                    done_rows.append(index-1)
                else:
                    continue
            elif index not in done_rows:
                last_list.append(row.tolist())
                done_rows.append(index)
        del_haplo_count=pd.DataFrame(last_list,columns=['Start','End','0','0.5','1','SNPS','Haplotype'])
        deleting_rows=[]
        # In order to join haploblocks after deleting small haploblocks, iterate over rows, comparing the row with the previous one.
        for index, row in del_haplo_count.iterrows():
            if index!=0:
                if del_haplo_count.loc[index]['Haplotype']==del_haplo_count.loc[index-1]['Haplotype']:  #If haploblocks have same phase
                    # Form an haploblock from both haploblocks
                    del_haplo_count.at[index,'Start']=del_haplo_count.loc[index-1]['Start']
                    del_haplo_count.at[index,'SNPS']=del_haplo_count.loc[index-1]['SNPS']+del_haplo_count.loc[index]['SNPS']
                    deleting_rows.append(index-1)
        del_haplo_count.drop(deleting_rows, inplace=True)
        del_haplo_count=del_haplo_count.reset_index() 
        del_haplo_count.drop('index',axis=1,inplace=True)
        return(del_haplo_count[['Haplotype', 'Start', 'End', 'SNPS']])


    def mark_snps(self,HetMatSign,i,cell):
        """ Function that groups raw consecutive SNPs with the same phase to form haploblocks."""
        list_phase=[]
        list_init_phase=[]
        list_end_phase=[]
        chromosome=[]
        previous_phase=float('nan')
        for index, row in HetMatSign[['Chr','Position',str(i)+'.'+str(cell)+'.Phase']].iterrows():
            if math.isnan(row[str(i)+'.'+str(cell)+'.Phase'])==False:   # Only use defined SNPs
                if math.isnan(previous_phase): # Detect if it is first SNP of the first haploblock
                    previous_phase=row[str(i)+'.'+str(cell)+'.Phase']
                    list_phase.append(row[str(i)+'.'+str(cell)+'.Phase'])
                    list_init_phase.append(row['Position'])
                    chromosome.append(row['Chr'])
                if row[str(i)+'.'+str(cell)+'.Phase']!=previous_phase:  # Detect if it is the first SNP of the haploblock
                    previous_phase=row[str(i)+'.'+str(cell)+'.Phase']
                    list_phase.append(row[str(i)+'.'+str(cell)+'.Phase'])
                    list_init_phase.append(row['Position'])
                    cond=True
                    previous_step=-1
                    while cond: # Find the previous SNP with a defined and different phase.
                        if math.isnan(HetMatSign.iloc[index+previous_step][str(i)+'.'+str(cell)+'.Phase']):
                            previous_step=previous_step-1
                        else:
                            list_end_phase.append(HetMatSign.iloc[index+previous_step]['Position'])
                            cond=False
            if index==HetMatSign[['Chr','Position',str(i)+'.'+str(cell)+'.Phase']].index.values[-1]: # Find the last defined SNP of the last haploblock
                cond=True
                previous_step=0
                while cond:
                    if math.isnan(HetMatSign.iloc[index+previous_step][str(i)+'.'+str(cell)+'.Phase']):
                        previous_step=previous_step-1
                    else:
                        list_end_phase.append(HetMatSign.iloc[index+previous_step]['Position'])
                        cond=False
        array_=np.array((list_phase,list_init_phase,list_end_phase))    #Return an array with the phase, the start and the end of the different haplotypes
        return(array_)

    def save_bed(self, dict, i):
        """ Function for saving the data in simple BED file per cell with chr, start coordinate, end coordinate and the phase"""
        if os.path.exists(str(self.output_dir)):
            pass
        else:
            os.makedirs(str(self.output_dir))
        for d in dict:
            dict[d]['Chr']= 'chr' + dict[d]['Chr'].astype(str)
            saving_file=pd.concat([dict[d]['Chr'], dict[d]['Start'],dict[d]['End'],dict[d]['Haplotype']], axis=1, keys=['Chr', 'Start','End','Haplotype'])
            saving_file.to_csv(str(self.output_dir)+'/Trio'+str(i)+'Cell'+str(d)+'.bed',index=False,header=False,sep='\t')

    def phasing(self):
        """Main function for phasing. Select only informative maternal SNPs. Choose egg in the first trio as reference.
        Do not consider heterozygous SNP in the reference or NC. Mark the phase of the SNPs. Apply the previous functions"""

        inputfile=self.inputfile
        # Select informative maternal SNPs
        HetMatIndex=inputfile.index[(inputfile['gDNA.GType']=='AB') | (inputfile['gDNA.GType']=='BA')]
        HetMatIndexPost=HetMatIndex+1
        HetMat=inputfile.iloc[HetMatIndex.tolist(),:]
        list_strings_egg=("\n".join(s for s in list(HetMat) if '.egg.GType'.lower() in s.lower()))
        list_num_trios=[int(s) for s in list(list_strings_egg) if s.isdigit()]
        # Ignore heterozygous or NC SNPs in the reference.
        HetMatSign=HetMat[(HetMat['1.egg.GType']!='NC')&(HetMat['1.egg.GType']!='AB')]
        HetMatSign=HetMatSign.reset_index()
        HetMatSign=HetMatSign.rename(columns = {'index':'SNP.num'})
        self.reference=HetMatSign['1.egg.GType']
        # Go through trios
        for i in list_num_trios:
            # Mark the phase of the SNPs
            # PB1 cell
            HetMatSign[str(i)+'.PB1.Phase']=pd.Series(float('nan'))
            HetMatSign.loc[:,str(i)+'.PB1.Phase'][(HetMatSign[str(i)+'.PB1.GType']!=self.reference)&(HetMatSign[str(i)+'.PB1.GType']!='NC')]=0. # Not in phase
            HetMatSign[str(i)+'.PB1.Phase'][(HetMatSign[str(i)+'.PB1.GType']==self.reference)&(HetMatSign[str(i)+'.PB1.GType']!='NC')]=1.   # In phase
            HetMatSign[str(i)+'.PB1.Phase'][(HetMatSign[str(i)+'.PB1.GType']=='AB')&(HetMatSign[str(i)+'.PB1.GType']!='NC')]=0.5    #Heterozygous
            # PB2 cell
            HetMatSign[str(i)+'.PB2.Phase']=pd.Series(float('nan'))
            HetMatSign[str(i)+'.PB2.Phase'][(HetMatSign[str(i)+'.PB2.GType']!=self.reference)&(HetMatSign[str(i)+'.PB2.GType']!='NC')]=0.   # Not in phase
            HetMatSign[str(i)+'.PB2.Phase'][(HetMatSign[str(i)+'.PB2.GType']==self.reference)&(HetMatSign[str(i)+'.PB2.GType']!='NC')]=1.   # In phase
            HetMatSign[str(i)+'.PB2.Phase'][(HetMatSign[str(i)+'.PB2.GType']=='AB')&(HetMatSign[str(i)+'.PB2.GType']!='NC')]=float('nan')   # Heterozygous: incorrect
            # EGG cell
            HetMatSign[str(i)+'.egg.Phase']=pd.Series(float('nan'))
            HetMatSign[str(i)+'.egg.Phase'][(HetMatSign[str(i)+'.egg.GType']!=self.reference)&(HetMatSign[str(i)+'.egg.GType']!='NC')]=0.   # Not in phase
            HetMatSign[str(i)+'.egg.Phase'][(HetMatSign[str(i)+'.egg.GType']==self.reference)&(HetMatSign[str(i)+'.egg.GType']!='NC')]=1.   # In phase
            HetMatSign[str(i)+'.egg.Phase'][(HetMatSign[str(i)+'.egg.GType']=='AB')&(HetMatSign[str(i)+'.egg.GType']!='NC')]=float('nan')   # Heterozygous: incorrect
            cells_results={}
            # Proceed with PB2 and EGG. The procedure of creating the haploblocks for these two types of cells is the same
            for cell in ['PB2','egg']:
                chr_results=[]
                # Iterate over Chromosomes
                for chr in HetMatSign['Chr'].drop_duplicates().values:
                    chr_HetMatSign=HetMatSign[HetMatSign['Chr']==chr]
                    array_=self.mark_snps(chr_HetMatSign,i,cell)    # Create simple (primitive) haploblocks
                    min_results=self.limit_snps(array_,chr_HetMatSign[[str(i)+'.'+str(cell)+'.Phase','Position']])  #Delete small haploblocks
                    final_results=min_results.sort_values(by=['Start']) #Final results for the chromosome
                    final_results['Chr']=chr
                    chr_results.append(final_results)
                cells_results[cell]=pd.concat(chr_results)  # Join all the chromosomes data for each cell
            chr_results=[]
            # Proceed with PB1. The procedure is a bit different due to the noise contained in the heterozygous haploblocks
            # Iterate over chromosomes
            for chr in HetMatSign['Chr'].drop_duplicates().values:
                chr_HetMatSign=HetMatSign[HetMatSign['Chr']==1]
                min_results=self.treat_PB1(cells_results,chr_HetMatSign,i)  # Apply the function to form the haploblocks in a cell where heterozygosity can happen 
                final_results=min_results.sort_values(by=['Start'])
                final_results['Chr']=chr
                chr_results.append(final_results)
            cells_results['PB1']=pd.concat(chr_results) #Join all the chrosomosomes data for each cell
            self.save_bed(cells_results,i)  # Save results for each trio
