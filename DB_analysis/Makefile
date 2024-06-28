.PHONY: run help plot rank cdb model new_db pred_ligands file_prep analysis 

help:
	@echo ''
	@echo -e '######################## HELP ##########################'
	@echo ''
	@echo -e '>>> Usage: diffduck [\033[1m\033[31moptions\033[0m]'
	@echo ''
	@echo -e '          \033[1m\033[31mrun\033[0m: Provided that you have placed the analyzed' 
	@echo -e '\tresults of the DiffDock run (whether each sample'
	@echo -e '\tattached to chain 1 or 0, i.e. results.csv) in'
	@echo -e '\tthe DD/ directory, generates all the outputs and'
	@echo -e '\tthe plots'
	@echo ''
	@echo -e '       \033[1m\033[31mnew_db\033[0m: Prints the column names of the database.'
	@echo ''
	@echo -e '          \033[1m\033[31mcdb\033[0m: Cleans the database and generates the necessary'
	@echo -e '\tfiles for generating a DiffDock input.'
	@echo ''
	@echo -e '         \033[1m\033[31mrank\033[0m: Ranks the results after the DiffDock run.'
	@echo ''
	@echo -e '         \033[1m\033[31mplot\033[0m: Plots the results after ranking the results.'
	@echo ''
	@echo -e '    \033[1m\033[31mfile_prep\033[0m: Generates the files and PDB necessary to'
	@echo -e '\trun the docking program.'
	@echo ''
	@echo -e '      \033[1m\033[31mdocking\033[0m: Uses the docking program selected to'
	@echo -e '\tpredict the ligand posses for the protein complex given.'
	@echo ''
	@echo -e '     \033[1m\033[31manalysis\033[0m: Analyzes the results of the docking program.'
	@echo ''
	@echo -e '        \033[1m\033[31mmodel\033[0m: Provides datafiles for trainning and generates'
	@echo -e '\ta model for ranking the DiffDock results. In order to'
	@echo -e '\tsolely run the models, data_A.csv and data_B.csv are'
	@echo -e '\tneeded. Use diffduck cdb in order to generate them.'
	@echo ''
	@echo -e '########################################################'
	@echo ''

cdb:
	@./codes/run_files/cdb.sh

rank:
	@./codes/run_files/rank.sh

plot:
	@./codes/run_files/plot.sh

model:
	@./codes/run_files/model.sh

new_db:
	@./codes/run_files/new_db.sh

file_prep:
	@./codes/run_files/file_prep.sh

docking:
	@./codes/run_files/sbatch.sh

analysis:
	@./codes/run_files/analysis.sh
	

run:cdb rank plot file_prep pred_ligands analysis model

