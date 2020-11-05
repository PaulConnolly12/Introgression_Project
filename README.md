A repository for scripting associated with analyzing differential divergence of cosmopolital Drosophila melanogaster populations to the ancestral Zambian population. UNDER ACTIVE DEVELOPMENT
---------------------------------------------------------------------------------------------
"Call me Ishmael. Some years ago—never mind how long precisely—having little or no money in my purse, and nothing particular to interest me on shore, I thought I would sail about a little and see the watery part of the world. It is a way I have of driving off the spleen and regulating the circulation. Whenever I find myself growing grim about the mouth; whenever it is a damp, drizzly November in my soul; whenever I find myself involuntarily pausing before coffin warehouses, and bringing up the rear of every funeral I meet; and especially whenever my hypos get such an upper hand of me, that it requires a strong moral principle to prevent me from deliberately stepping into the street, and methodically knocking people’s hats off—then, I account it high time to get to sea as soon as I can. This is my substitute for pistol and ball. With a philosophical flourish Cato throws himself upon his sword; I quietly take to the ship. There is nothing surprising in this. If they but knew it, almost all men in their degree, some time or other, cherish very nearly the same feelings towards the ocean with me."

Contributors
----------------------------------
Emma Howell
Paul Connoly
Chris McAllester


Available Scripts
----------------------------------

Divergence calculation: calculates the divergence in fixed or variable length windows across a set of population files organized as per-site tabbed counts of bases as 'A T G C N'
Script can be run from python, using the test data, as:

	python Scripts\divergence\divergence.py --d_mel_count_files testdata_CO_Chr2R.txt testdata_FR_Chr2R.txt --d_sim_ref_count_file testdata_sim_Chr2R.txt --d_mel_count_file_path test_data\ --d_sim_count_file_path test_data\ --window_size 10 --window_type total --out test_outputs\test_run 



Installation
----------------------------------

Download the scripts from this repository:

	git clone https://github.com/MesserLab/SLiM.git

Scripts present require python or R to run. 
