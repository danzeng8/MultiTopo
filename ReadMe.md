# Topological Simplification of Nested Shapes

This repository implements the paper found here: [https://github.com/danzeng8/MultiTopo/blob/main/nested_self.pdf](https://github.com/danzeng8/MultiTopo/blob/main/nested_self.pdf)

Experimental data can be found here: [https://drive.google.com/drive/folders/1QZRpfTymGDKu6qqkpyTQ5-1UN_uSAFOV?usp=sharing](https://drive.google.com/drive/folders/1QZRpfTymGDKu6qqkpyTQ5-1UN_uSAFOV?usp=sharing)

# Parameters:

--in : input file, which can be a sequence of image slices or a .tif file

--out : output file: which can also be image slices or a .tif file

--S : isovalues of the input image which are to be repaired. See below for examples. For time-series data, each time point is assigned to a particular iso-value in the input image

--indices : the indices of the shape (from --S) values which are to be repaired. See below for examples.

--inFileType : 0 (default) for image slices, 1 for .tif. If --in is a .tif file, change this to 1.

--distMode : determines how cuts and fills are computed. 0 is using intensity field (default), 1 is using uniform inflation / deflation (for time-series shapes)

--beam : width of beam search (default = 1)

--epsilon : parameter for epsilon simplification, used in the bone example. All cuts and fills with a maximum intensity difference greater than epsilon are not considered in the optimization (default = 0)

--propagate : Use the greedy propagation strategy instead of state-space search (default = false)

# Commands to reproduce paper results:

Figure 2: MultiTopo --S "0 1 2 3" --indices "1 2 3 4" --in examples_in/p0511/ --out examples_outp0511/ --distMode 1 

Figure 10: MultiTopo --in examples_in/angiography.tif --out examples_out/angiography_out.tif  --S "105 120 135 150 165 180" --indices "1 2 3 4 5 6" --inFileType 1

Figure 11: MultiTopo --S "40 60 80 100 120 140 160”  --indices "1 2 3 4 5 6 7" --in examples_in/bone/ --out examples_out/bone_out/ --epsilon 40 --distMode 1

Figure 12: MultiTopo --S "80 168 254" --in examples_in/brain/ --out examples_out/brain/ --distMode 1 --indices “1 2 3”

Figure 13: MultiTopo --S "0 1 2" --in examples_in/rice_root_3.tif --out examples_out/rice_root_out_3.tif --distMode 1

Figure 14: MultiTopo --S "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40" --indices "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41" --in examples_in/maize_root_42/  --out examples_out/maize_root_42/ --distMode 1

