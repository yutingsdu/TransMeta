Description
================

TransMeta is an assembler that can simultaneously assemble RNA-seq reads of multiple samples. TransMeta can output a unified set of meta-annotations(GTF format) for all samples in an RNA-seq experiment and use the information from all the samples to generate a set of transcripts for each individual sample. TransMeta constructs the novel vector-weighted splicing graph model, in which the edges and nodes are weighted by vectors rather than single numbers with the element in kth position of the vectors being the corresponding weight in sample k. And TransMeta employs a cosine-similarity based combing strategy and label setting algorithm to recover the transcripts. TransMeta takes a file that lists the alignment files(Bam format) of multiple RNA-seq samples as input, and outputs all assembled candidate transcripts in GTF format. 

The software has been developed to be user-friendly, which could be downloaded from https://github.com/yutingsdu/TransMeta

And it is free to use, modify, redistribute without any restrictions, except including the license provided with the distribution.


Installation
================

Prerequisites

 g++ with support for C++11 (e.g. 4.7.2)

 Boost >= 1.47.0

 Bamtools (now tested with 2.5.2 -- the makefile also contains instructions for older versions)

 libz

# 1. Installing Boost

    a) download a version of boost and unpack it

        $ tar zxvf boost_1_47_0.tar.gz

    b) change to the boost directory and run ./bootstrap.sh

        $ cd boost_1_47_0
        $ ./bootstrap.sh

    c) run

        $ ./b2 install --prefix=/your/boost/dir

        ########################################################################
        #  
        #  For example, 
        #
        #  if you want to install boost in /home/yuting/local/boost, the command is :
        #
        #  $ ./b2 install --prefix=/home/yuting/local/boost
        #
        #  If the boost is installed successfully, you would find two sub-directories:
        #
        #  /home/yuting/local/boost/include/
        #  /home/yuting/local/boost/lib/
        #
        #########################################################################

    Note: The default Boost installation directory is /usr/local. Take note of the boost 
    installation directory, because you need to tell the TransMeta installer where to find 
    boost later on.

# 2. Installing BamTools

    Download bamtools via: git clone https://github.com/pezmaster31/bamtools.git 

    Build bamtools by following the steps below.

    a) go to the bamtools directory and make a new directory named "build"

        $ mkdir build
        $ cd build

    b) type cmake and make it install

        $ cmake -DCMAKE_INSTALL_PREFIX=/your/bamtools/dir ..
        $ make
        $ make install

        where CMAKE_INSTALL_PREFIX is the root of your final installation directory.

        ##########################################################################
        #  
        # For example,
        # 
        # if you want to install bamtools in /home/yuting/local/bamtools, the command is :
        # 
        # $ cmake -DCMAKE_INSTALL_PREFIX=/home/yuting/local/bamtools ..
        # $ make
        # $ make install
        #    
        # If the bamtools is installed successfully, you would find two sub-directories:
        #
        # /home/yuting/local/bamtools/include
        # /home/yuting/local/bamtools/lib64
        #
        ##########################################################################

    As an altanitive, you can build bamtools based on the instruction at
    https://github.com/pezmaster31/bamtools/wiki/Building-and-installing
    
    Note: you need to tell the TransMeta installer where to find bamtools later on.

# 3. Building TransMeta

    Change to the TransMeta-v.1.0/src directory and make

        $ cd src
        $ make all BOOST_PATH=/your/boost/dir BAMTOOLS_PATH=/your/bamtools/dir

        where BOOST_PATH is the aformentioned directory where you installed the boost 
        and BAMTOOLS_PATH is the directory where you installed the bamtools.

        #########################################################################
        #
        # For example,
        #
        # if you installed boost in /home/yuting/local/boost and 
        # installe bamtools in /home/yuting/local/bamtools, the command is :
        #
        # $ make all BOOST_PATH=/home/yuting/local/boost BAMTOOLS_PATH=/home/yuting/local/bamtools
        #
        ##########################################################################

    If the TransMeta is installed successfully, you'll see the following 7 executable files in TransMeta-v.1.0/src/bin/
    transmeta_abundance, transmeta_cover, transmeta_graph, transmeta_merge, transmeta_path_search, transmeta_individual, and transmeta_AG_assembly.

# 4. Running TransMeta

    a) Type the following command OR Set the LD_LIBRARY_PATH enviroment variable

        $ export LD_LIBRARY_PATH=/home/yuting/local/boost/lib:$LD_LIBRARY_PATH

        Note: please replace "/home/yuting/local/boost/lib" with your own directory "/your/boost/dir/lib"

        ##########################################################################
        #           !!!!!!!!!! PLEASE NOTE !!!!!!!!!!
        #
        #  If you do not set this variable , you would possible see the follwoing error information:
        #
        # "error while loading shared libraries: libboost_serialization.so.1.47.0:
        #                  cannot open shared object file: No such file or directory"
        #
        ##########################################################################

    b) The executable TransMeta is in the TransMeta-v.1.0 directory

        $ TransMeta -B bamFile_list -s first -o TransMeta_outdir -p 10

        where the bamFile_list is a file that lists the alignments BAM files (one per line) of multiple samples.


# 5. Testing TransMeta on a demo data set:

    To test if you have succesfully installed TransMeta, please download the demo data set from 
      
    https://sourceforge.net/projects/transmeta/files/DemoData/. 

    At this website you will see two alignments files of two samples(sample1.bam and sample2.bam)

    Put the Hisat.bam and Star.bam in the directory TransMeta-v.1.0/sample_test/ and change to TransMeta-v.1.0/sample_test/

    Type the following command:

        $ ./run_me.sh

    If you get the TransMeta_outdir/TransMeta.gtf, congratulations, you have succesfully installed the TransMeta.



===========================================================================

    TransMeta v.1.0 usage:

    ** Required **

    --bam/-B <string>		: path to the file listing the alignments BAM files (one per line)

    --strand/-s <string> 		: Strand-specific RNA-Seq reads orientation.

			   If reads are paired:
				    1) Use <unstranded> to indicate RNA-seq reads are non-strand-specific.
				    2) Use <first> to indicate fr-first-stranded RNA-seq reads.
				    3) Use <second> to indicate fr-second-stranded RNA-seq reads.

			   If reads are single:
				    1) Use <single_unstranded> to indicate RNA-seq reads are non-strand-specific.
				    2) Use <single_forward> to indicate RNA-seq reads are forward.
				    3) Use <single_reverse> to indicate RNA-seq reads are reverse.

---------------------------------------------------------------------------

    ** Options **

    --help/-h			  : Output TransMeta Help Information.

    --version/-v			  : Print current version of TransMeta.

    --output_dir/-o <string>	  : Output path, default: transmeta_outdir.

    --min_meta_trans_cov/-c <float>	  : Minimum expression level estimated by abundance analysis for the meta_assembly, default: 2.0.
				  (Please Note: TransMeta automatically outputs assemblies with different expression levels in the output_dir)

    --min_indiv_trans_cov/-i <float>  : Minimum expression level estimated by abundance analysis for the certain assembly of each sample, default: 1.5.

    --annotation_reference/-g <string>: reference annotation to use for guiding the assembly process (GTF).

    --min_trans_length/-L <int>   	  : Minimum assembled transcript length, default: 500.

    --min_average_frac/-d <float> 	  : Minimum junction coverage fraction by average junction coverage, default: 0.03.

    --min_unbalance_frac/-D <float>   : Minimum fraction of unbalanced junction, default: 0.03.

    --thread/-p <int> 		  : Number of threads to use default: 1.

---------------------------------------------------------------------------

** Typical commands **

  (i) A typical TransMeta command for paired-end data might be:

    TransMeta -B bamFiles_list -s first -o TransMeta_outdir -p 25

  (ii) A typical TransMeta command for single-end data might be:

    TransMeta -B bamFiles_list -s single_reverse -o TransMeta_outdir -p 25

===========================================================================


Authors: Ting Yu designed and wrote TransMeta.

Contact:
 Any questions, problems, bugs are welcome and should be dumped to
 Ting Yu <yutingsdu@163.com>
 Created on Nov 24, 2021.

