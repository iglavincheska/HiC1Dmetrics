import pandas as pd
import numpy as np
import os
import random
import subprocess
import sys

def loadDenseMatrix(filename,log=False):
    #print(filename)
    if not os.path.exists(filename):
        print("File does not exist, please check your file path")
        exit(1)

    print("loading data...")
    try:
        data = pd.read_csv(filename, delimiter='\t', index_col=0)
    except:
        print('Error input matrix, please check https://h1d.readthedocs.io/en/latest/overview.html#input-format')
        exit(1)
    print("loading finished")

    if log == True:
        return(np.log1p(data))
    else:
        return data

def loadWithNorm(filename,method= "RPM",log = False):
    print("loading data...")
    data = pd.read_csv(filename, delimiter='\t', index_col=0)
    if method == "RPM":
        data = (10000000 * data) / np.nansum(data)
    print("loading finished")

    if log:
        return np.log1p(data)
    else:
        return data

def hic2matrix(path,res,chr,gt,juicer=None):
    if not gt:
        print("rawhic require genome_table file");exit(1)
    try:
        gtfile = pd.read_csv(gt,sep="\t",header=None)
        if not int(gtfile.iloc[0,1]):
            exit(1)
    except:
        print("Wrong genome_table file.")
        print("Please check your genome_table file. Is it tab separated ?")
        exit(1)

    print("Start dump matrix from hic file")
    codepath = os.path.dirname(os.path.realpath(__file__))
    makeIntra = codepath+"/extract/makeMatrixIntra.sh"

    if not juicer:
        print("Use built-in juicertool v1.11.04")
        juicer = codepath+"/jc/jctool_1.11.04.jar"
    else:
        print("Use customized juicertool")

    foldername = "./MatrixTemp"+str(random.random())
    os.system("bash "+makeIntra+" "+"KR"+" "+"."+" "+path+" "+str(res)+" "+gt+" "+juicer+" "+chr+" "+foldername + "> info.txt")
    matrixpath = foldername+"/"+str(res)+"/observed.KR."+chr+".matrix.gz"
    print("Finish dump")
    return(matrixpath)

def cool2matrix(path,res,chr,gt):
    print("Start dump matrix from cool file")
    if not gt:
        print("cool require genome_table file");exit(1)
    codepath = os.path.dirname(os.path.realpath(__file__))
    makeIntra = codepath+"/extract/coolerdump.sh"
    foldername = "./MatrixTemp"+str(random.random())
    cmd = "H1D_PYTHON='{}' bash '{}' '{}' '{}' '{}' '{}' '{}' > info.txt".format(
        sys.executable, makeIntra, path, str(res), chr, gt, foldername
    )
    ret = os.system(cmd)
    matrixpath = foldername+"/"+str(res)+"/"+chr+".matrix.gz"
    if ret != 0 or (not os.path.exists(matrixpath)):
        print("Error: cool2matrix failed. Please check info.txt and your Python environment.")
        exit(1)
    print("Finish dump")
    return(matrixpath)

def cool2hic(path, res, gt, juicer=None):
    """
    Convert .cool or .mcool file to .hic format using juicer_tools.
    
    Args:
        path: Path to .cool or .mcool file
        res: Resolution for the output .hic file
        gt: Path to genome table file
        juicer: Path to juicer_tools jar file (optional, uses built-in if not specified)
    
    Returns:
        Path to the generated .hic file
    """
    if not gt:
        print("cool2hic requires genome_table file"); exit(1)
    
    try:
        import cooler
    except ImportError:
        print("Error: cooler package is required for .cool/.mcool conversion")
        print("Install with: pip install cooler")
        exit(1)
    
    # Validate genome table
    try:
        gtfile = pd.read_csv(gt, sep="\t", header=None)
        if not int(gtfile.iloc[0, 1]):
            exit(1)
    except:
        print("Wrong genome_table file.")
        print("Please check your genome_table file. Is it tab separated?")
        exit(1)
    
    print("Start converting cool file to hic format")
    
    # Detect file type
    file_type = path.split('.')[-1]
    if file_type == 'cool':
        c = cooler.Cooler(path)
    elif file_type == 'mcool':
        c = cooler.Cooler(path + '::resolutions/' + str(res))
    else:
        print("Error: Input file must be .cool or .mcool"); exit(1)
    
    # Extract chromosomes (exclude those with "_" or "M")
    chrom = c.chroms()[0:len(c.chroms())]
    chrom = chrom.sort_values('name')
    
    # Create temporary matrix file in juicer format
    temp_matrix = "./cool_temp_" + str(random.random()) + ".txt"
    
    if os.path.exists(temp_matrix):
        os.remove(temp_matrix)
    
    sel = c.matrix(balance=False, as_pixels=True, join=True, sparse=True)
    chroms = list(chrom.name)
    
    print("Converting cool matrix to juicer format...")
    for a_idx, chrA in enumerate(chroms):
        for chrB in chroms[a_idx:]:
            # chrA==chrB gives cis; chrA!=chrB gives trans
            hic_df = sel.fetch(chrA, chrB)
            
            if hic_df.empty:
                continue
            
            # keep only needed columns
            hic_df = hic_df.loc[:, ["chrom1", "start1", "chrom2", "start2", "count"]].copy()
            
            # add dummy juicer fields
            hic_df["str1"] = 0
            hic_df["frag1"] = 0
            hic_df["str2"] = 0
            hic_df["frag2"] = 1
            
            # order for juicer short-with-score style
            hic_df = hic_df[["str1", "chrom1", "start1", "frag1",
                             "str2", "chrom2", "start2", "frag2", "count"]]
            
            hic_df.to_csv(temp_matrix, sep="\t", mode="a", index=False, header=False)
    
    print("Juicer format matrix created at:", temp_matrix)
    
    # Find juicer tools
    codepath = os.path.dirname(os.path.realpath(__file__))
    if not juicer:
        juicer = codepath + "/jc/jctool_1.11.04.jar"
        if not os.path.exists(juicer):
            juicer = codepath + "/jc/juicer_tools.2.20.00.jar"
        print("Use built-in juicertool")
    else:
        print("Use customized juicertool")
    
    if not os.path.exists(juicer):
        print("Error: Juicer tools jar not found at", juicer); exit(1)
    
    # Convert matrix to .hic using juicer pre
    hic_output = path.replace("." + file_type, "") + "_KR.hic"
    
    print("Converting to .hic format using juicer_tools...")
    cmd = f"java -Xms512m -Xmx20480m -jar {juicer} pre -r {res} -k KR {temp_matrix} {hic_output} {gt}"
    print("Running:", cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print("Error converting to .hic:", result.stderr)
        exit(1)
    
    # Clean up temporary matrix file
    if os.path.exists(temp_matrix):
        os.remove(temp_matrix)
    
    print("Finish converting to .hic")
    return hic_output
#===========================
