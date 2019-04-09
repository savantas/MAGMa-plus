MAGMa+
======

* MAGMa+ allows MAGMa metabolite identification with dynamic parameter selection.
* Disclaimer: The installation of MAGMa/MAGMa+ was only tested on a linux machine.
* Prerequisites: A working version of Python2.7 is required
* Remark: The MAGMa version provided here was originally written by: Lars Ridder

Installation of MAGMa+
----------------------

* Make sure you have all MAGMa dependencies installed
 * the most important dependency is RDKit (http://www.rdkit.org/) with INCHI support.
 * other required python libaries are: sqlalchemy, lxml, numpy, pp, request, macauthlib, mock, nose, coverage, sklearn

* Copy all the MAGMa/MAGMa+ source files and accompanying files to a dedicated directory.


Running MAGMa+
--------------

* Setting up the shell environment: `MAGMAPLUS_CLASSIFIER_PATH`

	If the classifier files and AllFingerprints.txt file are in a different directory than the one you run MAGMa+ from,
	the `MAGMAPLUS_CLASSIFIER_PATH` variable should be set to that different directory.
	It should contain: `AllFingerprints.txt`, `metabolite_posD_rfc.pkl` and `metabolite_negC_rfc.pkl`

* Unzip the classifier pickle files (`metabolite_posD_rfc.pkl.zip` and `metabolite_negC_rfc.pkl.zip`).

* Creating a configuration file

  A `magma_job.ini` config file is read from working directory.

  Example content:

```
	[magma job]
	# Location of structure database to fetch candidate molecules to match against ms peak trees
	# db is expected to be available at where job is executed
	structure_database.hmdb = /path/to/HMDB_MAGMa.db
	chemical_engine = rdkit
```

* Generating a molecular structure database.

  A script (`process_hmdb.py`) is provided that generates an HMDB database. It can be adapted at will to generate other structure databases.

* Run MAGMa+ exactly as you would run MAGMa. For more information, type:

```
	python /path_to_magma_plus/MAGMa_plus.py -h
```

* Example:

```
	export MAGMAPLUS_CLASSIFIER_PATH=/path/to/classifiers
	echo '100.112725211: 100 (53.0388: 0.836306820171, 55.0544: 25.3576989727, 58.029: 0.108793463482, 83.0855: 100.0)' > example.tree
	python MAGMa_plus.py read_ms_data -i 1 -p 5 -q 0.001 -f mass_tree example.tree output.db
	python MAGMa_plus.py annotate -c 0 -d 0 -b 3 -w 1 -s hmdb output.db
	python MAGMa_plus.py export_result output.db
```

Citing MAGMa+
-------------
Please cite the papers:
- http://link.springer.com/article/10.1007%2Fs11306-016-1036-3
- http://onlinelibrary.wiley.com/doi/10.1002/rcm.6364/full
