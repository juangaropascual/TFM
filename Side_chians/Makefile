.PHONY: help clean ddi sbatch_docking analysis structures plot mutate structures_sbatch rank_all

help:
	@echo ''
	@echo -e '######################## HELP ##########################'
	@echo ''
	@echo -e '>>> Usage: peligand [\033[1m\033[31moptions\033[0m]'
	@echo ''
	@echo -e '    \033[1m\033[31mfold\033[0m: generates the structures with AlphaFold, DeppFold ' 
	@echo -e '\t& OmegaFold.'
	@echo ''
	@echo -e '     \033[1m\033[31mddi\033[0m: generates an input for the performance of the docking.' 
	@echo ''
	@echo -e ' \033[1m\033[31mdocking\033[0m: Uses the docking program selected to predict the'
	@echo -e '\tligand posses for the protein complex given.'
	@echo ''
	@echo -e '\033[1m\033[31manalysis\033[0m: Analyzes the results of the docking program.'
	@echo ''
	@echo -e '    \033[1m\033[31mrank\033[0m: ranks the results of the docking simulations.' 
	@echo ''
	@echo -e '    \033[1m\033[31mplot\033[0m: Plots the ranked results of the simulations.' 
	@echo ''
	@echo -e '  \033[1m\033[31mmutate\033[0m: Performs n_mutations single aa mutations to the protein.' 
	@echo ''
	@echo -e '########################################################'
	@echo ''

docking:
	@./codes/run_files/docking_sbatch.sh

analysis:
	@./codes/run_files/analysis_sbatch.sh

ddi:
	@./codes/run_files/c_i.sh

fold:
	@./codes/run_files/structures_sbatch.sh

rank:
	@python codes/rank_all.py

mutate:
	@./codes/run_files/mutate.sh

plot:
	@./codes/run_files/plot.sh

search:
	@python codes/search_DB.py	

clean: 
	@rm -rf protein_structures/*
	@rm -rf results/*
	@rm -f inputs/docking_*
	@mkdir results/figures
