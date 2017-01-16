# CC4.0
Cross-compare script with rewritten algorithm. Able to execute in more than one process.

### What is CC?
CC stands for cross-compare. It's a python script originally written to empirically evaluate datasets for optimally orthogonal submatricies. Data is read in from ``.csv`` files and can be processed on multiple threads using python's ``multiprocessing`` module. For large datasets or submatrix dimensions (>2) it is advisable to use a supercomputing cluster to run this algorithm.

### Dependencies
The following must be installed in order for the code to run properly:  

- Python3 ([https://www.python.org/]())
- Pandas ([http://pandas.pydata.org/]())
- Numpy ([http://www.numpy.org/]())

For easy installation and many more useful python modules try the [Anaconda installer](https://www.continuum.io/downloads).


### Running the script

**Formatting the input:**  
Data must be formatted as a ``.csv`` file with row 1 as column headings, and column 1 as row indices. Each entry in the table must be a number, or may be blank. Strings will break the code. Place substrates in the columns of the  A small table might look like this:  

|			|Substrate 1	|Substrate 2	|Substrate 3	|...	|Substrate n
|---		|---			|---			|---	|---	|---
|Enzyme 1	|3.4			|5.7			|10	|...	|...
|Enzyme 2	|67				|8.9			|11	|...	|...
|Enzyme 3	|4.8			|10E10			|5	|...	|...
|...		|...			|...			|...	|...	|...
|Enzyme m	|...			|...			|...	|...	|value for substrate n in enzyme m

**Running your ``.csv``:**  
After forking, cloning, or downloading the repository, navigate to the folder containing the source files using the terminal.

On mac/linux: ``cd path/to/CC4.0``

A basic run can be achieved by specifying only the input and output filenames with the ``-i`` and ``-o`` tags respectively.

On mac/linux: ``python3 run_OSF.py -i "name of your.csv" -o "name of your desired output.csv"``

The script will print what it is doing into the terminal and create a ``.csv`` in the current path with the specified file name.

**Other options for customizing the script output:**  
``[-d DIMENSION]`` Optional tag to specify the number of dimensions to use in the search. The default is 2. Using greater than 2 dimensions on large (>500 entry) datasets could take a considerable amount of time.

``[-p PROCESSES]`` Number of child processes to split the iteration into. If a personal computer is being used, putting the number of available processors here will give the speediest result. The default is 1.

``[-l LENGTH]`` Length of the list to be outputted. Putting ``-l 2000`` will give the top 2000 ranked submatricies. Default: 1000.

``[-t THRESHOLD]`` RMS threshold to keep for the sorted list. Specifying this tag is often necessary for larger searches (>2 dimensions or larger datasets) to reduce memory usage. Setting this to 0.15 reduced our memory usage adequately. Default: 1.

**Example input:**  
For a 3 dimensional run with a file called ``dataset.csv``, on a computer with 8 cores, one may use:

``python3 run_OSF.py -i "dataset.csv" -o "dataset_out.csv" -d 3 -t 0.1 -p 8``
