# Change this as required.
# Note that the path shouldn't contain spaces! (GNU Make is old!!)
# eventually the destination path should be passed as env so no folder renaming is necessary
SOLUTION_PATH=/home/launchpad/out/variables/bin/sol.
OUTDATA_PATH=/home/launchpad/out/variables/bin
GNU_CSPLIT=csplit
GCC=gcc

csvs: bin_outdata/bin_csvs/xc.csv
csv_archive: compressed_csvs.zip

.PHONY: clean
clean:
	rm -rf split python_arrays parsed_solution.txt export_solution

compressed_csvs.zip: bin_outdata/bin_csvs/xc.csv
	zip $@ $(OUTDATA_PATH)/*/*.csv

bin_outdata/bin_csvs/xc.csv: split/file00 parse_variables.py
	python3 parse_variables.py ./split/ $(OUTDATA_PATH)

split/file00: parsed_solution.txt
	mkdir -p split
	$(GNU_CSPLIT) --quiet -z -f split/file $^ /Contents*/ '{*}'

parsed_solution.txt: export_solution $(SOLUTION_PATH)mds
	./export_solution $(SOLUTION_PATH) | tee $@

export_solution: export_solution.c
	$(GCC) -o $@ export_solution.c

export_solution_debug: export_solution.c
	$(GCC) -o $@ export_solution.c -g
