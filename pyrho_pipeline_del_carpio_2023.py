# Workflow for the pyrho analyses
# Created by Maria Izabel Cavassim Alves and modified by Tina Del Carpio

# This workflow requires the following software to be installed in your conda environment:
# bcftools
# smc++
# pyrho

#this workflow is a pipeline for gwf manual: https://gwf.app/
#for a guide on using gwf on UCLA's HPC (Hoffman2) see: 
#https://github.com/izabelcavassim/gwf_UCLA_cluster/blob/master/gwf_tutorial.md

from gwf import Workflow
import glob
import os
import os.path
import itertools
#import pandas as pd

gwf = Workflow()


def merge_vcf_files(directory, species, size, outdir):
	inputs = []
	outputs=[f'{outdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_allchr.vcf.gz']
	options = {
	'memory':'8g',
	'cores':8,
	'walltime':'23:59:59',
	}
	spec = f'''
	source /u/local/Modules/default/init/modules.sh
	module load bcftools

	bcftools concat /u/scratch/c/cad17/wolf_recomb/dog_genomes_clares_cassdata_final_filtered_vcf_justSNPs_6files/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr{{1..39}}.vcf.gz -Oz -o {outdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_allchr.vcf.gz
	'''
	# print(spec)
	# print(outputs)
	return inputs, outputs, options, spec

#this function will conver your vcf files to smc format for smc++

def demography_file_prep(species, size, inputdir, outdir, roh):
	inputs = [f'{inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr{i}.vcf.gz' for i in list(range(1, 38))]
	if roh == False:
		outputs=[f'{outdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr{i}.smc.gz' for i in list(range(1, 38))]
		#print(outputs)
	else:
		outputs=[f'{outdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_roh_masked_chr{i}.smc.gz' for i in list(range(1, 38))]
	options = {
	'memory':'8g',
	'cores':8,
	'walltime':'23:59:59',
	#smc++ vcf2smc {input_vcf} {inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_allchr.smc.gz

	}
	input_vcf = inputs[0]
	if roh == False:
		spec = f'''
		/u/local/Modules/default/init/modules.sh
		. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
		conda activate /u/home/m/mica20/.conda/envs/noncoding

		for i in {{1..39}};
			do smc++ vcf2smc {inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr$i.vcf.gz {outdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr$i.smc.gz chr$i {species}:{species}1,{species}2,{species}3,{species}4,{species}5,{species}6,{species}7,{species}8,{species}9,{species}10,{species}11,{species}12,{species}13,{species}14,{species}15,{species}16,{species}17,{species}18 --ignore-missing
		done

		'''
	else:
		spec = f'''
		/u/local/Modules/default/init/modules.sh
		. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
		conda activate /u/home/m/mica20/.conda/envs/noncoding

		for i in {{1..39}};
			do smc++ vcf2smc {inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr$i.vcf.gz {outdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_roh_masked_chr$i.smc.gz -c 2000000 chr$i {species}:{species}1,{species}2,{species}3,{species}4,{species}5,{species}6,{species}7,{species}8,{species}9,{species}10,{species}11,{species}12,{species}13,{species}14,{species}15,{species}16,{species}17,{species}18,{species}19,{species}20,{species}21,{species}22,{species}23,{species}24,{species}25  --ignore-missing
		done

		'''
	if species == "AW" and roh == True:
		spec = f'''
		/u/local/Modules/default/init/modules.sh
		. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
		conda activate /u/home/m/mica20/.conda/envs/noncoding

		for i in {{1..39}};
			do smc++ vcf2smc {inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr$i.vcf.gz {outdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_roh_masked_chr$i.smc.gz -c 2000000 chr$i {species}:{species}10,{species}11,{species}12,{species}13,{species}14,{species}15,{species}16,{species}17,{species}18,{species}19,{species}20,{species}21,{species}22,{species}23,{species}24,{species}25  --ignore-missing
		done

		'''
	if species == "AW" and roh == False:
		spec = f'''
		/u/local/Modules/default/init/modules.sh
		. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
		conda activate /u/home/m/mica20/.conda/envs/noncoding

		for i in {{1..39}};
			do smc++ vcf2smc {inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr$i.vcf.gz {outdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr$i.smc.gz -c 2000000 chr$i {species}:{species}10,{species}11,{species}12,{species}13,{species}14,{species}15,{species}16,{species}17,{species}18,{species}19,{species}20,{species}21,{species}22,{species}23,{species}24,{species}25  --ignore-missing
		done

		'''
	#print(spec)
	return inputs, outputs, options, spec

#this function infers the demography with smc++
def demography_estimation(species, inputdir, size, roh=True):	
	if roh == False:
		inputs=[f'{inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr{i}.smc.gz' for i in list(range(1, 38))]

		outputs=[f'/u/home/m/mica20/project-kirk-bigdata/canid_project/scripts/output_{species}/model.final.json', f'/u/home/m/mica20/project-kirk-bigdata/canid_project/scripts/plot_{species}.pdf',f'/u/home/m/mica20/project-kirk-bigdata/canid_project/scripts/plo\
t_{species}.csv' ]

	else:
		inputs=[f'{inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_roh_masked_chr{i}.smc.gz' for i in list(range(1, 38))]

		outputs=[f'/u/home/m/mica20/project-kirk-bigdata/canid_project/scripts/output_{species}/model.final.json', f'/u/home/m/mica20/project-kirk-bigdata/canid_project/scripts/plot_{species}_roh.pdf',f'/u/home/m/mica20/project-kirk-bigdata/canid_project/scripts/plo\
t_{species}_roh.csv' ]
	options = {
	'memory':'8g',
	'cores':8,
	'walltime':'23:59:59',
	}
	if roh == False:
		spec = f'''
		/u/local/Modules/default/init/modules.sh
		. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
		conda activate /u/home/m/mica20/.conda/envs/noncoding

		smc++ estimate 4.5e-9 -o output_{species}/ {inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr*
		
		smc++ plot plot_{species}.pdf output_{species}/model.final.json -c 
		'''
	if roh == True:
		spec = f'''
		/u/local/Modules/default/init/modules.sh
		. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
		conda activate /u/home/m/mica20/.conda/envs/noncoding

		smc++ estimate 4.5e-9 -o output_{species}/ {inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_roh_masked_chr*
		
		smc++ plot plot_{species}_roh.pdf output_{species}/model.final.json -c 
		'''

	#print(spec)
	return inputs, outputs, options, spec

#this function was used to split our vcf files by chromosome, but can be skipped if your files already are
def spliting_vcf_files(species, dir_files, chrom, sample):
	''' Spliting VCF files per chromosome to run with pyrho '''
	inputs = [f'{dir_files}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{sample}{species}_jointcalled_allchr.vcf.gz']
	outputs=[f'{dir_files}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{sample}{species}_jointcalled_chr{chrom}.vcf.gz', f'{dir_files}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{sample}{species}_jointcalled_chr{chrom}.vcf.gz.tbi']
	
	options = {
	'memory':'8g',
	'cores':8,
	'walltime':'23:59:59',
	}
	bcftools view -S {individuals_file_directory}{species}_indivs.txt {inputs[0]} --force-samples | bgzip -c > {dir_files}/chroms/{species}/{sample}{species}_chr{chrom}.vcf.gz
	vcftools  --gzvcf  {inputs[0]}  --chr chr{chrom}  --recode --recode-INFO-all --out  {dir_files}/chroms/{species}/{sample}{species}_chr{chrom}.vcf
	bcftools view -r 1 file.vcf.gz

	individuals_file_directory = '/u/home/m/mica20/project-kirk-bigdata/canid_project/aw_pg_vcfs_and_demogs/'
	spec = f'''
	/u/local/Modules/default/init/modules.sh
	. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
	conda activate /u/home/m/mica20/.conda/envs/noncoding

	bcftools view -r chr{chrom} {inputs[0]} | bgzip -c > {dir_files}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{sample}{species}_jointcalled_chr{chrom}.vcf.gz
	tabix -p vcf {dir_files}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{sample}{species}_jointcalled_chr{chrom}.vcf.gz
	'''
	#print(spec)
	return inputs, outputs, options, spec

#this part defines functions for pyrho specifically

#function to generate a lookup table of values of rho and allele frequencies
def make_lookup_table(n, N, mu, smcpp_file, species, chrom, results_dir):
	''' Lookuptable for pyrho ''' 
	#inputs = [f'{dir_files}/Results/{smcpp_file}']
	inputs = [f'{smcpp_file}']
	outputs=[f'{results_dir}/logfile_{species}_n_{n}_N_{N}_lookuptable_{chrom}_{day}.txt', f'{results_dir}/{chrom}_{species}_n_{n}_N_{N}_lookuptable_{day}.hdf']
	options = {
	'memory':'50G',
	'cores':8,
	'walltime':'23:59:59',
	'user':'cad17',
	'email':'bea',
	}
	spec = f'''
	source /u/local/Modules/default/init/modules.sh
	module load gcc/4.9.3
	module load anaconda3
	source $CONDA_DIR/etc/profile.d/conda.sh
	conda activate pyrho
	module load gsl
	module load hdf5
	
	pyrho make_table -n {n} -N {N} --mu {mu} --logfile {results_dir}/logfile_{species}_n_{n}_N_{N}_lookuptable_{chrom}_{day}.txt --outfile {results_dir}/{chrom}_{species}_n_{n}_N_{N}_lookuptable_{day}.hdf \
	--approx --smcpp_file {smcpp_file} --decimate_rel_tol 0.1
	'''.format(n=n, N=N, mu=mu, smcpp_file=smcpp_file, chrom=chrom, species=species, results_dir=results_dir)
	print(spec)
	return inputs, outputs, options, spec

#function that tests window size and block penalty
def get_hyperparam(n, N, mu, smcpp_file, species, chrom, sims, dir_files, results_dir):
	'''  Getting the hyperparam ''' 
	inputs = [f'{results_dir}/{chrom}_{species}_n_{n}_N_{N}_lookuptable_{day}.hdf', f'{smcpp_file}'] 
	outputs=[f'{results_dir}/{chrom}_{species}_n_{n}_N_{N}_hyperparam_results_{day}.txt']
	options = {
	'memory':'8g',
	'cores':1,
	'walltime':'23:59:59',
	'user':'cad17',
	'email':'bea',
	}
	spec = f'''
	source /u/local/Modules/default/init/modules.sh
	module load gcc/4.9.3
	module load anaconda3
	source $CONDA_DIR/etc/profile.d/conda.sh
	conda activate pyrho
	module load gsl
	module load hdf5

	pyrho hyperparam -n {n} --mu {mu} --blockpenalty 50,100 \
	--windowsize 25,50 --logfile {results_dir}/logfile_{species}_n_{n}_N_{N}_hyperparameter_{day}.txt  --tablefile {results_dir}/{chrom}_{species}_n_{n}_N_{N}_lookuptable_{day}.hdf \
	--num_sims {sims} \
	--smcpp_file {smcpp_file} --outfile {results_dir}/{chrom}_{species}_n_{n}_N_{N}_hyperparam_results_{day}.txt 
	'''.format(n=n, N=N, mu=mu, smcpp_file=smcpp_file, chrom=chrom, sims=sims, results_dir=results_dir)
	print(spec)
	return inputs, outputs, options, spec

#function to estimate recombination with pyrho
def optimize(n, N, mu, species, chrom, vcf, blockpenalty, windowsize, results_dir, dir_files):
	'''  Optmizing parameters (blockpenalty and windowsize) ''' 
	inputs = [f'{results_dir}/{chrom}_{species}_n_{n}_N_{N}_lookuptable_{day}.hdf', f'{vcf}']
	outputs=[f'{results_dir}/logfile_{species}_n_{n}_N_{N}_optimize_{chrom}_{day}.txt', f'{results_dir}/{chrom}_{species}_n_{n}_N_{N}_{day}.rmap']
	options = {
	'memory':'8g',
	'cores':1,
	'walltime':'23:59:59',
	'user':'cad17',
	'email':'bea',
	}
	spec = f'''
	source /u/local/Modules/default/init/modules.sh
	module load gcc/4.9.3
	module load anaconda3
	source $CONDA_DIR/etc/profile.d/conda.sh
	conda activate pyrho
	module load gsl
	module load hdf5
	
	pyrho optimize --tablefile {results_dir}/{chrom}_{species}_n_{n}_N_{N}_lookuptable_{day}.hdf \
	--logfile {results_dir}/logfile_{species}_n_{n}_N_{N}_optimize_{chrom}_{day}.txt \
	--vcffile {vcf} \
	--ploidy 2 \
	--outfile {results_dir}/{chrom}_{species}_n_{n}_N_{N}_{day}.rmap \
	--blockpenalty {blockpenalty} --windowsize {windowsize} \
	'''.format(n=n, N=N, mu=mu, chrom=chrom, blockpenalty=blockpenalty, windowsize=windowsize, results_dir=results_dir)
	print(spec)
	return inputs, outputs, options, spec


#this section was used to take the specific formatting of our file names and make 
#easier to read in systematically


Making a dictionary of breeds and sample sizes
metadata = {"AW":15, "BC":7, "EW":10, "IR": 10, "LB":10, "PG":15, "TM":9} # "AW" is 15, but is actually 14
metadata = {"BC":"Ten", "EW":"Ten", "IR":"Ten", "LB":"Ten", "TM":"Ten", "AW":"Fifteen", "PG":"Fifteen"} #, "MD":"group", "MW":"group"} # "AW" is 15, but is actually 14

dir_files = "/u/home/m/mica20/project-kirk-bigdata/canid_project/aw_pg_vcfs_and_demogs"
results_dir = "/u/scratch/c/cad17/wolf_recomb/Results"

dir_files = "/u/scratch/c/cad17/wolf_recomb"
chroms = range(1,39,1)
species = ["AW", "PG"]
metadata = {"AW":15, "BC":7, "EW":10, "IR": 10, "LB":10, "PG":15, "TM":9} # "AW" is 15, but is actually 14
metadata = {"BC":"Ten", "EW":"Ten", "IR":"Ten", "LB":"Ten", "TM":"Ten", "AW":"Fifteen", "PG":"Fifteen"} #, "MD":"group", "MW":"group"} # "AW" is 15, but is actually 14

day = "2023_10_01"

#from here on is where we fed the details of our files to actually run the functions defined above 

#####################################################
# Subsetting files                                  #
#####################################################
dir_files = "/u/scratch/c/cad17/wolf_recomb/dog_genomes_clares_cassdata_final_filtered_vcf_justSNPs_6files"
for sp in species:
	for ch in chroms:
		gwf.target_from_template(f'Subsetting_vcf_{sp}_{ch}', spliting_vcf_files(dir_files=dir_files, chrom=ch, species=sp, sample=metadata[sp]))


######################################################
# COMPUTING DEMOGRAPHY.                              #
######################################################
smcpp = True
if smcpp == True:
	species = metadata.keys()
	vcf_dir =  "/u/scratch/c/cad17/wolf_recomb/dog_genomes_clares_cassdata_final_filtered_vcf_justSNPs_6files"
	out_dir = "/u/scratch/c/cad17/wolf_recomb/chroms/CONCATENATED"
	for sp in species:
		print(sp)
		vcf_dir_sp = vcf_dir
		size = metadata[sp]

		spliting_vcf_files
		gwf.target_from_template(f'Making_smcpp_input_{sp}', demography_file_prep(species=sp, size=size, inputdir=vcf_dir_sp, outdir=out_dir, roh=True))
		gwf.target_from_template(f'Estimating_demography_{sp}',demography_estimation(species=sp, size=size, inputdir=out_dir, roh=True))

		###
		## Estimating stats for each vcf file
		gwf.target_from_template(f'Merging_vcf_{sp}', merge_vcf_files(directory=vcf_dir_sp, species=sp, size=size, outdir=out_dir))
		gwf.target_from_template(f'Producing_VCF_stats_{sp}', vcf_stats(sample=size, species=sp))
		#gwf.target_from_template(f'Producing_VCF_stats_PG', vcf_stats(dir_files=dir_files, sample=15, species='PG'))


######################################################
# COMPUTING Recombination  Pyrho                     #
######################################################
pyrho = True
metadata = {"PG":15}
metadata_names = {"PG":"Fifteen"}
species = metadata.keys()
roh = False
chroms = range(1,39,1)


if pyrho == True:
	#for loop here
	for sp in species:
		for ch in chroms:
			#gwf.target_from_template(f'Subsetting_vcf_{sp}_{ch}', spliting_vcf_files(dir_files=dir_files, chrom=ch, species=sp, sample=metadata[sp])) #this is if you need to split your vcfs by chromosome

			# Making lookup tables for each chromosome
			lower_string = sp.lower()
			number_string = metadata_names[sp]
			vcf_name = f'BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{number_string}{sp}_jointcalled_chr{ch}.vcf.gz'
			smcpp_dir = "/u/home/c/cad17/project-klohmuel/wolf_recomb/smc_files" #edit this
			gwf.target_from_template(f'Producing_lookuptable_{sp}_{ch}', make_lookup_table(n=2*metadata[sp], N=round(2*metadata[sp]+ (metadata[sp]*2)*0.5), mu=4.5e-9, smcpp_file=f'{smcpp_dir}/{sp}_snps_all_autosomes_demog_ibdne_wolfNa_2023_09_22.csv', species=sp, chrom=ch, results_dir=results_dir))	
			gwf.target_from_template(f'getting_hyperparam_{sp}_{ch}', get_hyperparam(n=2*metadata[sp], N=round(2*metadata[sp]+ (metadata[sp]*2)*0.5), mu=4.5e-9, sims=10, smcpp_file=f'{smcpp_dir}/{sp}_snps_all_autosomes_demog_ibdne_wolfNa_2023_09_22.csv', species=sp, chrom=ch, dir_files=dir_files, results_dir=results_dir))
			sample = metadata[sp]
			gwf.target_from_template(f'optimize_parameters_{sp}_{ch}', optimize(n=2*metadata[sp], N=round(2*metadata[sp]+ (metadata[sp]*2)*0.5), mu=4.5e-9, vcf=f'{dir_files}/{sp}/{vcf_name}', species=sp, chrom=ch, blockpenalty=50, windowsize=50, results_dir=results_dir, dir_files=dir_files))
			if roh == True:
				gwf.target_from_template(f'Producing_lookuptable_{sp}_{ch}_roh', make_lookup_table(n=2*metadata[sp], N=round(2*metadata[sp]+ (metadata[sp]*2)*0.5), mu=4.5e-9, smcpp_file=f'{smcpp_dir}/plot_{sp}_roh.csv', species=sp, chrom=ch, dir_files=dir_files, results_dir=results_dir))	
				gwf.target_from_template(f'getting_hyperparam_{sp}_{ch}_roh', get_hyperparam(n=2*metadata[sp], N=round(2*metadata[sp]+ (metadata[sp]*2)*0.5), mu=4.5e-9, sims=10, smcpp_file=f'{smcpp_dir}/plot_{sp}_roh.csv', species=sp, chrom=ch, dir_files=dir_files, results_dir=results_dir))
				sample = metadata[sp]
				gwf.target_from_template(f'optimize_parameters_{sp}_{ch}_roh', optimize(n=2*metadata[sp], N=round(2*metadata[sp]+ (metadata[sp]*2)*0.5), mu=4.5e-9, vcf=f'{dir_files}/{sp}/{vcf_name}', species=sp, chrom=ch, blockpenalty=50, windowsize=50, results_dir=results_dir, dir_files=dir_files))
