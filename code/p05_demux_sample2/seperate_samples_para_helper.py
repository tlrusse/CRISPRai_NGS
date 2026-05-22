#! /home/craut/miniconda3/envs/mageck/bin/python
# this python script is used to seperate the samples into the different barcodes
# this script has the following command line arguments:
# -s or --samfile: the sam file to be analyzed
# -c or --cigar: a comma-separated string of cigar strings to be kept from the sam file (default="6S30M56S")
# -o or --output: the output file name (will be stored as output.fastq)
# -n or --ncores: the number of cores to use for parallel processing (default=number of cores on the machine -1)
# -d or --debug: whether or not to print debug statements (default=False)
# -b or --buffer_size: the buffer size in lines to read from the sam file at a time (default=1_000_000)
# example: ./seperate_samples_para_helper.py -s test.sam -c "6S30M115S" -o ./demuxed/test_sam_dm -n 4 -d

import sys, os
import shutil
from pathlib import Path
import argparse
from multiprocessing import Pool, cpu_count
from typing import Set, Tuple, List, Optional
from tqdm import tqdm
from contextlib import ExitStack

# set and grab the cli arguments
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--samfile", help="the sam file to be analyzed", required=True)
parser.add_argument("-c", "--cigar", help="a comma-separated string of cigar strings to be kept from the sam file (default='6S30M56S')", default="6S30M56S")
parser.add_argument("-o", "--output", help="the output file name (will be stored as output_RNAME.fastq)", required=True)
parser.add_argument("-n", "--ncores", help=f"the number of cores to use for parallel processing (default={cpu_count()-1})", type=int, default=cpu_count()-1)
parser.add_argument("-d", "--debug", help="whether or not to print debug statements (default=False)", action="store_true")
parser.add_argument("-b","--buffer_size", help="the buffer size in lines to read from the sam file at a time (default=1_000_000)", type=int, default=1_000_000)
args = parser.parse_args()

IO_BUFFER_SIZE_FOR_SHUTIL = 16 * 1024 * 1024  # 16 MB buffer size for shutil.copyfileobj, which is a good balance between speed and memory usage

def calculate_read_offsets(sam_file:str, ncores:int, debug:bool = False) -> list:
    """_summary_

    Args:
        sam_file (str): the sam file to be analyzed
        ncores (int): the number of cores to use for parallel processing
        debug (bool, optional): whether or not to print debug statements. Defaults to False.

    Returns:
        list: a list of int tuples containing the start and end byte offsets for each core's read in the sam file
    """
    # calculate the total size of the sam file
    total_size = os.path.getsize(sam_file)
    if debug:
        print(f"total_size: {total_size}")
    # handle the edge case when ncores is 1 or 0 or greater than the number of bytes in the sam file
    if ncores == 0:
        raise ValueError("Error: ncores cannot be 0")
    elif ncores == 1:
        return [(0, total_size)]
    elif ncores > total_size:
        raise ValueError("Error: ncores cannot be greater than the number of bytes in the sam file. Is the sam file empty?")
    
    # calculate the size of each chunk for each core    chunk_size = total_size / ncores
    chunk_size = total_size // ncores
    if debug:
        print(f"chunk_size: {chunk_size}")
    
    # calculate the byte offsets for each core's read in the sam file
    offsets = []
    for i in range(ncores):
        start = i * chunk_size
        end = start + chunk_size if i < ncores - 1 else total_size
        offsets.append((start, end))
    if debug:
        print(f"initial offsets: {offsets}")
    # check that each tuple in offsets if within +/- ncores of the expected chunk size
    # and that the start and end offsets are the same across 2 adjacent tuples
    for i in range(ncores - 1):
        if abs(offsets[i][1] - offsets[i][0] - chunk_size) > ncores:
            raise ValueError(f"Warning: chunk size for core {i} is not within +/- {ncores} of the expected chunk size")
        if offsets[i][1] != offsets[i + 1][0]:
            raise ValueError(f"Warning: end offset for core {i} does not match start offset for core {i + 1}")
    # check last chunk
    if abs(offsets[-1][1] - offsets[-1][0] - chunk_size) > ncores:
        raise ValueError(f"Warning: chunk size for core {ncores - 1} is not within +/- {ncores} of the expected chunk size")
    return offsets

def create_intermediate_directories(output_header: str, ncores: int, debug: bool = False) -> None:
    """_summary_

    Args:
        output_header (str): the output file name header (will be stored as output_RNAME.fastq)
        ncores (int): the number of cores to use for parallel processing
        debug (bool, optional): whether or not to print debug statements. Defaults to False.
    """
    if ncores <= 0:
        raise ValueError("Error: ncores must be a positive integer greater than 0.")

    # Initialize the base intermediate directory path object
    base_dir = Path(f"{output_header}_parallel_helpers")
    
    # Check existence beforehand ONLY for logging purposes if debug is enabled
    if debug:
        if base_dir.exists():
            print(f"Intermediate directory already exists: {base_dir}")
        else:
            print(f"Created intermediate directory: {base_dir}")
            
    # parents=True creates missing parent folders; exist_ok=True prevents crashes if it exists
    base_dir.mkdir(parents=True, exist_ok=True)

    # Generate isolated subdirectories for each core
    for i in range(ncores):
        core_dir = base_dir / f"core_{i}"  # The '/' operator joins paths beautifully in pathlib
        
        if debug:
            if core_dir.exists():
                print(f"Intermediate directory for core {i} already exists: {core_dir}")
            else:
                print(f"Created intermediate directory for core {i}: {core_dir}")
                
        core_dir.mkdir(parents=True, exist_ok=True)

def delete_intermediate_directories(output_header: str, debug: bool = False) -> None:
    """Delete the intermediate directories created for parallel processing.

    Args:
        output_header (str): the output file name header
        debug (bool, optional): whether or not to print debug statements. Defaults to False.
    """
    base_dir = Path(f"{output_header}_parallel_helpers")
    
    if base_dir.exists():
        # shutil.rmtree deletes everything inside, completely safely and recursively
        shutil.rmtree(base_dir)
        
        if debug:
            print(f"Deleted intermediate directory: {base_dir}")
    else:
        if debug:
            print(f"Intermediate directory does not exist, so nothing to delete: {base_dir}")

def construct_parallel_args(sam_file: str, output_header: str, cigars: Set[str], offsets: list, buffer_size: int, debug: bool = False) -> list:
    """Construct the list of arguments to be passed to the parallel processing function.

    Args:
        sam_file (str): the sam file to be analyzed
        output_header (str): the output file name header (will be stored as output_RNAME.fastq)
        cigars (set): a set of cigar strings to be kept from the sam file
        offsets (list): a list of int tuples containing the start and end byte offsets for each core's read in the sam file
        buffer_size (int): the buffer size in lines to read from the sam file at a time
        debug (bool, optional): whether or not to print debug statements. Defaults to False.

    Returns:
        list: a list of tuples containing the arguments for each core's processing function
    """
    args = []
    base_dir = Path(f"{output_header}_parallel_helpers")
    # extract the fname starter from the output header. Get by subtracting the directory from the output header
    fname_starter = os.path.basename(output_header)
    for i, (start, end) in enumerate(offsets):
        # The filenames stay the same across the cores so merging downstream is easier. 
        # Results are segregated by the intermediate directories.
        worker_output_header = base_dir / f"core_{i}" / fname_starter
        args.append((sam_file, str(worker_output_header), cigars, start, end, i, buffer_size))
        if debug:
            print(f"Constructed arguments for core {i}: {(sam_file, worker_output_header, cigars, start, end, i, buffer_size)}")
    return args

def process_line(line: str, cigars: Set[str], output_header: str) -> Optional[Tuple[str, str]]:
    """Function to process a single line from the sam file 

    Args:
        line (str): a single line from the sam file
        cigars (set): a set of cigar strings to be kept from the sam file
        output_header (str): the output file name header (will be stored as output_RNAME.fastq)
    """
    if line.startswith("@"):
        # this is a header line, so we can skip it
        return None
    fields = line.strip().split("\t",12)
    if len(fields) < 11:
        return (line.strip()+ "\n", f"{output_header}_incomplete.sam")
    qname, rname, cigar, seq, qual = fields[0], fields[2], fields[5], fields[9], fields[10].rstrip()
    
    ofname = f"{output_header}_{rname}.fastq" if rname != "*" and cigar in cigars else f"{output_header}_unmapped.fastq"
    
    return (f"@{qname}\n{seq}\n+\n{qual}\n", ofname)

def empty_buffer(buffer: list) -> None:
    """Function to write the contents of the buffer to the appropriate files and clear the buffer.

    Args:
        buffer (list): a list of tuples containing the lines to be written and the corresponding file names
    """
    # first we want to group the lines by the file name
    file_dict = {}
    for line, file_name in buffer:
        if file_name not in file_dict:
            file_dict[file_name] = []
        file_dict[file_name].append(line)
    # now we can write the lines to the appropriate files
    for file_name, lines in file_dict.items():
        with open(file_name, "a") as f:
            f.write("".join(lines))
    buffer.clear()

def process_chunk(args: Tuple[str, str, Set[str], int, int, int, int]) -> None:
    """Function to process a chunk from the sam file and demux the reads into separate fastq files based on the RNAME and CIGAR string.
    
    Args:
        args (tuple): a tuple containing the following arguments:
            sam_file (str): the sam file to be analyzed
            output_header (str): the output file name header
            cigars (set): a set of cigar strings to be kept from the sam file
            start_bytes (int): the starting byte offset for this chunk
            end_bytes (int): the ending byte offset for this chunk
            core_id (int): the id of the core processing this chunk (used for naming intermediate files)
            buffer_size (int): the buffer size in lines to read from the sam file at a time
    """

    sam_file, output_header, cigars, start_bytes, end_bytes, core_id, buffer_size = args
    # open the sam file and seek to the starting byte offset
    with open(sam_file, "r") as f:
        f.seek(start_bytes)
        # if this is not the first chunk, we need to read until we reach a newline character to ensure we start at the beginning of a line
        if core_id != 0:
            _ = f.readline() # discard the first line since it may be incomplete
        # now we can read the lines until we reach the ending byte offset
        buffer = []
        while f.tell() < end_bytes:
            line = f.readline()
            processed_line = process_line(line, cigars, output_header)
            if processed_line:
                buffer.append(processed_line)
            # if the buffer size is reached, write the buffer to the appropriate files and clear the buffer
            if len(buffer) >= buffer_size:
                empty_buffer(buffer)
        # after the loop, we may have some remaining lines in the buffer that need to be written to the appropriate files
        if buffer:
            empty_buffer(buffer)

def merge_intermediate_files(folders: List[str], output_folder: str, debug: bool = False) -> None:
    """Merge the fastq and sam files from the specified folders into the output_folder, ensuring that files with the same name are concatenated together.

    Args:
        folders (List[str]): a list of folders containing the intermediate files to be merged
        output_folder (str): the folder where the merged files will be stored
        debug (bool, optional): whether or not to print debug statements. Defaults to False.
    """
    fnames_in_folders = [set(os.listdir(folder)) for folder in folders]
    all_fnames = set().union(*fnames_in_folders)
    for fname in all_fnames:
        nmerges = 0
        output_path = os.path.join(output_folder, fname)
        # Open a context stack for this specific barcode file
        with ExitStack() as stack:
            of = None
            for folder in folders:
                cur_fname = os.path.join(folder, fname)
                # check if the cur_fname exists
                if not os.path.exists(cur_fname):
                    continue
                if of is None:
                    # we want to just move the first file we encounter with this name to the output folder to avoid unnecessary copying
                    shutil.move(cur_fname, output_path)
                    # now we want to open the output file in append binary mode to write to it for the remaining files
                    of = stack.enter_context(open(output_path, "ab"))
                else:
                    # we want to append the contents of the current file to the output file
                    with open(cur_fname, "rb") as f:
                        shutil.copyfileobj(f, of, length=IO_BUFFER_SIZE_FOR_SHUTIL)
                nmerges += 1
        # exit stack automatically closes of 
        if debug:
            print(f"Merged file: {output_path}, Number of merges: {nmerges}")

def parallel_process_chunks(args_list: list, ncores: int) -> None:
    """Function to parallel process the chunks from the sam file using the multiprocessing Pool.

    Args:
        args_list (list): a list of tuples containing the arguments for each core's processing function
        ncores (int): the number of cores to use for parallel processing
    """

    with Pool(processes=ncores) as pool:
        list(tqdm(pool.imap_unordered(process_chunk, args_list), total=len(args_list), desc="Processing chunks in parallel"))

def main(s:str,c:Set[str],o:str,ncores:int,buffer_size:int,debug:bool = False) -> None:
    """_summary_

    Args:
        s (str): sam file to be analyzed
        c (set): set of cigar strings to be kept from the sam file
        o (str): output filename header
        ncores (int): number of cores to use for parallel processing
        buffer_size (int): the size of the buffer to use for copying files
        debug (bool, optional): whether or not to print debug statements. Defaults to False.
    """
    if debug:
        print(f"samfile: {s}")
        print(f"cigars: {c}")
        print(f"output_header: {o}")
        print(f"ncores: {ncores}")
        print(f"buffer_size: {buffer_size}")
    
    # calculate the byte offsets for each core's read in the sam file
    offsets = calculate_read_offsets(s, ncores, debug)
    # now we want to create the intermediate directories for each core
    create_intermediate_directories(o, ncores, debug)
    # construct the list of arguments to be passed to the parallel processing function
    args_list = construct_parallel_args(s, o, c, offsets, buffer_size, debug)
    # parallel process the chunks from the sam file using the multiprocessing Pool
    parallel_process_chunks(args_list, ncores)
    # after all the parallel processing is done, we want to merge the intermediate files into the output folder
    intermediate_folders = [str(Path(f"{o}_parallel_helpers/core_{i}")) for i in range(ncores)]
    # extract the output folder from the output header
    output_folder = os.path.dirname(o) if os.path.dirname(o) else "."
    merge_intermediate_files(intermediate_folders, output_folder, debug)
    # after merging, we can delete the intermediate directories
    delete_intermediate_directories(o, debug)

if __name__ == "__main__":
    sam_file = args.samfile
    cigars = set(args.cigar.split(","))
    output_header = args.output
    ncores = args.ncores
    buffer_size = args.buffer_size
    debug = args.debug
    main(sam_file,cigars,output_header,ncores,buffer_size,debug)
    