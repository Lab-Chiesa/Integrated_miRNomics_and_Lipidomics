# author: stefano.manzini@gmail.com

__version__ = "0.3.17.071024"
__author__ = "Stefano Manzini; stefano.manzini@gmail.com"

"""
> This library requires numpy, scipy and related scientific packages.
All of it is packaged into Anaconda - https://www.continuum.io/downloads

Also, it's better/required that you additionally install:

> Seaborn - https://stanford.edu/~mwaskom/software/seaborn/

pip install seaborn

> Plotly - https://plot.ly/python/getting-started/

pip install plotly

> Alternate regex module - https://pypi.python.org/pypi/regex

Download it, then cd into the dir and:
<your python> setup.py install

/* 07/10/24
 * Now strip_mirname() also works with the repo's pre-release files
 */

/* 23/06/23
 * well, well, plotly is deprecated. Gotta remove it.
 */
 

/* 05/11/2018
 * added make_data()
 */

/* 31/10/2018
 * added despace()
 */

/* 04/08/2018
 * added lipid_id_to_metaboanalyst(lipid_id)
 */

/* 20/06/2018
 * added absolute_filename()
 */

/* 16/05/2018
 * pre-compiled regexes for lipi_ids() and mirna_ids()
 */

/* 13/04/2018
 * added attrib_index()
 */

/* 23/11/2017
 * Reconciler-only functions moved over to Reconciler
 */

/* 31/10/2017
 * changed get_lipid_class2 to also evaluate free fatty acids of Cariplo 2.1
 */

/* 18/7/2017
 * added listify_GEO_genesymbols()
 */

/* 27/6/2017
 * added normalize_range()
 */

/* 23/5/2017
 * added uncomment()
 */

/* 28/9/2016:
 * - bugfix in:
 * mir_stdev_test()
 * lip_stdev_test()
 * When 0 averages && standard deviations are found in the two supplied
 * lists, now the test continues instead of returning True.
 *
 * - merged fuctions from forked library versions
 */

/* 27/9/2016
 * - bugfix in:
 * stdev_test()
 * When 0 averages && standard deviations are found in the two supplied
 * lists, now the test continues instead of returning True.
 */
"""

try:
    import regex as re
except:
    import re

import pandas as pd
import time
import os

from scipy.stats import (pearsonr, # (<pearson’s coefficient>, <2-tail-p-value)
                         spearmanr
                         )
                         
from numpy import (average,
                  polyfit,
                  poly1d,
                  linspace,
                  absolute,
                  sum)

from random import choice, randint, gauss

# 23/06/23
try:
    tls.set_credentials_file(username = "Manz", api_key = "2gmj2v0hxj")
except:
    pass


# =============================================================================
#                                * functions *
# =============================================================================

def make_data(average, stdev, n_datapoints, decimals=2):
    """Leverages on random.gauss to easily print datapoints to quick
    copy-and-paste operations
    """
    made_data = [gauss(average, stdev) for _ in range(n_datapoints)]
    for datapoint in made_data:
        print(round(datapoint, decimals))


def despace(stringlike):
    """Return a list of words separated by one or more space.
    It was difficult to make it with regex! :)
    Useful to humanize R's ANOVA Tukey post hocs and put them into spreadsheets
    """                                                
    returnlist = []
    word = ""
    for char in stringlike:
        if char == " ":
            if len(word) > 0:
                returnlist.append(word)
                word = ""
        else:
            word += char
    return returnlist


def lipid_id_to_metaboanalyst(lipid_id):
    """This function takes a Zora-produced lipid ID and
    returns a MetaboAnalyst-compatible lipid ID.

    This way, lipids are readable in MetaboAnalyst outputs.


    input: string
    returns: string
    """

    new = lipid_id.replace(
        ":","_"
        ).replace(
        "("," "
        ).replace(
        ")"," "
    )
    return new


def absolute_filename(fullpath):
    """This function strips a full path of all directories, leaving only the
    filename:
    
    /dir1/dir2/dir3/file

    becomes

    file

    input: string
    returns: string

    """
    match = re.search(r"(.*\/)(.*)", fullpath)
    return match.group(2)


def attrib_index(series):
    """Reads a column of a DataFrame, spits a list of indexes in order of
    appearance.
    
    Example considering a <DataFrame> df: 
    
                   TREND  
    gene                             
    Snap25         Lower in Western  
    Ptprz1         Lower in Western  
    6330403A02Rik  Higher in Western  
    Casq1          Higher in Western  
    Mylk4          No change  
    Acta1          No change  
    Slc6a2         No change  
    Ngfr           No change  
    Dbh            Higher in Bl/6  
    Myh2           Higher in Bl/6
    
    >>> attrib_index(df["TREND"])
    
    returns:
    
    [1, 1, 2, 2, 3, 3, 3, 3, 4, 4]
    
    Returns: <list>
    """
    returnlist = []
    keeptrack = []
    index = 0
    for elem in series:
        if elem not in keeptrack:
            index += 1
            keeptrack.append(elem)
            returnlist.append(index)
        else:
            returnlist.append(index)
    return returnlist


# ============================================================
def listify_GEO_genesymbols(dataframe, header = "Gene.symbol"):
# ============================================================
    """
    This function reads the "Gene.title" column out of a pandas dataframe
    generated by the online GEO top 250 analysis.

    This list is likely flawed by pandas NaN values, as well as genes identified as
    FAM69A///SNORA66///RPL5

    NaN values will be discarded, and genes resolved to individual names.

    This function returns a *new <list>*
    """
    returnlist = []
    try:
        for elem in dataframe[header]:
            try:
                for gene in elem.split("///"):
                    returnlist.append(gene)
            except:
                pass # (discards NaN)
    except KeyError:
        print("Can't find Gene.symbol: try specifying another column name via: header = \"new_header\"")
        return False
    except TypeError:
        print("This function accepts dataframes only. Trying to use it with lists?")
        return False
        
    return returnlist




# ========================
def normalize_range(
    iterable,
    allowed_max = 100,
    #allowed_min = 0,
    force_conversion = True
    ):
# ========================
    """Takes a range of numbers, and returns a modified list of items
    (same order) stretched or compressed to fit between 0  and defined max.
    
    force_conversion tries to floatize string-like numbers
    """
    
    if force_conversion is True:
        iterable = [float(x) for x in iterable]
    
    for n in iterable:
        if n < 0:
            raise ValueError("Not able to work with negative values, yet")
    
    mx = max(iterable)
    
    scale_factor = allowed_max / mx
    
    new_iterable = [x * scale_factor for x in iterable]
    
    return new_iterable
    
# // end of normalize_range //



# =========================================
def uncomment(stringlike, *, comment = "#"):
#?=========================================
    """Uncomments a Python-style commented line.
    Returns False if the entire line is commented,
    or the uncommented text otherwise.

    Examples:

    # hello there!
    returns False

    Ano1 # this is the Anoctamin gene
    returns "Ano1"

                  # comment starting after whitespace
    returns False
    """

    if stringlike.startswith("#"):
        return False
    elif len(stringlike.rstrip()) == 0:
        return False # empty line
    else:
        if "#" not in stringlike:
            return stringlike.rstrip().rstrip("\n")
        else:
            match1 = re.search(r"^.[^#]*", stringlike)
            returnstring = match1.group(0).rstrip().rstrip("\n")  # this strips any
                                                     # ending whitespace
            if len(returnstring) == 0:
                return False
            else:
                return returnstring
# // end of uncomment() //


# !Cariplo 2.1!
# ======================
def normalize_by_class(
        number,
        lipid_class,
        decimals = 0,
):
# ======================
    """Returns a number, divided by the number of class members.
    Supplying a lipid class is mandatory.
    """
    items_in_each_class = {'CE': 42, 'Cer': 17, 'DAG': 26, 'FC': 1, 'Gb3': 9,
                           'Glc': 10, 'Lac': 10, 'LPC': 19, 'LPE': 4, 'LPI': 1,
                           'PA': 3, 'PC': 75, 'PC O': 1, 'PC P': 8, 'PE': 52,
                           'PE O': 4, 'PE P': 20, 'PG': 10, 'PI': 7, 'PS': 18,
                           'SM': 26, 'TAG': 24}

    divisor = items_in_each_class[lipid_class]
    return str(round(number/divisor, decimals))


# ===================================================
def humanize_mirname(stringlike, cut_species = True):
# ===================================================
    """
    Takes something like:
    mmu-miR-434-5p_ID=MIM0001421_MI0001526_mmu-mir-434

    Spits something like (if cut_species is True):
    miR-434-5p

    Otherwise, it spits
    mmu-miR-434-5p
    """

    match = re.search(r"(.+?)_", stringlike)
    if cut_species is True:
        return match.group(1)[4:]
    else:
        return match.group(1)



# =================
def strip_extension(
    filename,
    extensions = (".txt", ".csv", ".tdt", ".xls", ".xlsx"),
    recursive = False,  # the first extension found, is stripped and that's it
    aggressive = False # tries to remove all extensions more aggressively
):
# =================
    
    """Strips a given filename <string> of the extension, if the filename
    ends in one such extensions.
    The first extension encountered is taken away, and the filename is
    returned. To strip all extension(s) recursively, set recursive to True.
    
    *please note* that the extension are removed in the order they appear
    in the given list. To try and remove all of them regardless of the order,
    set the function to aggressive mode, until I think a better way to
    implement this (I'm in a hurry now, and it's pointless to remove multiple
    extensions anyway).
    """
    
    if aggressive is False: 
        for extension in extensions:
            if filename.endswith(extension):
                if recursive is True:
                    filename = filename.rstrip(extension)
                else:
                    return filename.rstrip(extension)
    
        return filename
    else:   # aggressive = True
        for i in range((len(extensions) * len(extensions))):
            for extension in extensions:
                if filename.endswith(extension):
                    if recursive is True:
                        filename = filename.rstrip(extension)
                    else:
                        return filename.rstrip(extension)
    
        return filename



# ========================
def count_lines(filename):
# ========================
    
    """Counts the number of lines in a texftile named <filename>.
    <filename> is expected as a string, *not* as a <File object>.
    """
    
    counter = 0
    with open(filename, "r") as i:
        for line in i:
            counter += 1
    
    return counter


# =============================
def generate_random_key(length):
# =============================
    """Gives back a random text of given <int> length"""
    keymap = "abcdefghijklmnopqrstuvwxyz"
    returnstring = ""
    
    for i in range(length):
        char = choice(keymap)
        if randint(0, 100) > 50:
            returnstring += char
        else:
            returnstring += char.upper()
    
    return returnstring



# =======================================
def how_many_times_exceeds(df, threshold):
# =======================================
    
    """Horribly unoptimized function that counts the number of times
    a value exceeds the defined threshold in a pandas dataframe.
    """
    
    rows_number = df.shape[0]
    cols_number = df.shape[1]
    
    counter = 0    
    for row in range(rows_number):
        for col in df.columns:
            try:
                if df[col][row] >= threshold:
                    counter += 1
                else:
                    continue
            except:
                print("something went wrong: weird value (?):", df[col][row])
                raise
    
    return counter



# =======================================
def how_many_times_is_less(df, threshold):
# =======================================
    
    """Horribly unoptimized function that counts the number of times
    a value is lower than the defined threshold in a pandas dataframe.
    """
    
    rows_number = df.shape[0]
    cols_number = df.shape[1]
    
    counter = 0    
    for row in range(rows_number):
        for col in df.columns:
            try:
                if df[col][row] <= threshold:
                    counter += 1
                else:
                    continue
            except:
                print("something went wrong: weird value (?):", df[col][row])
                raise
    
    return counter




# ======
def ls(*,
       shorten_filenames = True,
       filenames_chunks = 20,
       pwd_at_the_end = True,
       file_id = "     ", # replaces <DIR> when printing files
       normalize_file_lenght = True,
       separator = "   "
):
# ======
    """ Kind of mimics ls bash command"""
    curr_path = os.getcwd()

    if shorten_filenames == True:
        dir_entries = [x[:filenames_chunks] + ".." + x[-filenames_chunks:] if len(x) > (filenames_chunks* 2) + 1 else x for x in os.listdir(".")]
    else:
        dir_entries = os.listdir(".")
    
    if normalize_file_lenght == True:
        dir_entries = [x if len(x) >= (filenames_chunks * 2) else x + (" " * ((filenames_chunks * 2) - len(x) + 2)) for x in dir_entries ]

    dir_entries_bool = [os.path.isdir(x) for x in os.listdir(".")]
    dir_entries_type = ["<DIR>" if x is True else file_id for x in dir_entries_bool]

    dir_entries_unixtime = [os.path.getmtime(x) for x in os.listdir(".")]
    dir_entries_time = [time.ctime(int(x)) for x in dir_entries_unixtime]

    dir_entries_byte_size = [os.path.getsize(x) for x in os.listdir(".")]
    total_file_size = sum(dir_entries_byte_size)

    if not pwd_at_the_end:
        print(f"\n Directory of {curr_path}")
    
    for file, ftype, ftime, fsize  in zip(
            dir_entries,
            dir_entries_type,
            dir_entries_time,
            dir_entries_byte_size
            ):
        print(file,
              separator,
              ftype,
              separator,
              ftime,
              separator,
              "{:,}".format(fsize) + " bytes",
              sep = ""
              )

    entries_number = len(dir_entries)
    dir_number = dir_entries_type.count("<DIR>")
    file_number = entries_number - dir_number

    if pwd_at_the_end:
        print(f"\nDirectory of {curr_path}")

    if entries_number > 1:
        print("\t\t", str(entries_number), "Entries")
    else:
        print("\t\t", str(entries_number), "Entries")

    if file_number > 1 or file_number == 0:
        print("\t\t", str(file_number), "Files",
              "\t{:,}".format(total_file_size), "bytes")
    else:
        print("\t\t", str(file_number), "File",
              "\t{:,}".format(total_file_size), "bytes")

    if dir_number > 1 or file_number == 0:
        print("\t\t", str(dir_number), "Directories")
    else:
        print("\t\t", str(dir_number), "Directory")



# ========================================
def index_lipid_species(dataframe, header):
# ========================================
    
    """Takes a Pandas <DataFrame> as input, and looks in specified
    header for lipids. It then builds a <dict>, with lipid classes as
    <keys> and lists of all lipid of that class as <values>.
    
    Returns: a <dict> (keys = lipid classes, values = all lipids in that class)
    """
    
    # get_lipid_class2()
    
    mm = {}
    
    for lipid in dataframe[header]:
        
        curr_class = get_lipid_class2(lipid)
        mm.setdefault(curr_class, [])
        mm[curr_class].append(lipid)
    
    return mm
        


# ====================================
def description_to_gene_symbol(string):
# ====================================
    
    
    """Takes an Ncbi description, spits back the gene symbol.
    Returns None if this cannot be accomplished.
    
    Example 1:
    
    description = "This gene is really (really) awesome (Aws1), mRNA."
    
    returned gene symbol = "Aws1"
    
    Example 2:
    
    description = "This example is much more common (Cmmn2), mRNA."
    
    returned gene symbol = "Cmmn2"
    """
    
    # looking for genes after one parenthesis (very very rare)
    try:
        
        match = re.search(r"\(.*\).*(\(.*\))", string)
        gene = match.group(1)
        return gene
    
    # looking for normally-looking gene descriptions
    except:
        
        try:
            
            match = re.search(r"(\(.*\))", string)
            gene = match.group(1)
            return gene
            
        except:
            
            return None



# =====================================================
def delta_percent_from_average(listlike, decimals = 18):
# =====================================================

    """ This function evaluates the % difference from each value of given
    list from the average of all values, and *returns* a list of such %.

    i.e.

    input = (10, 15, 20)    # the average is 15

    returned = [-0.50, 0, 0.25]
    """

    result = []
    av = average(listlike)
    for n in listlike:
        result.append(round(float((n - av)/n), decimals))

    return result



# ============
def pcorr_test(
        df,
        row = 2,
        threshold = 70
        ):
# ============
    
    """This returns Bool after evaluating,
    in a DataFrame (likely produced by Reconciler), if the Pearson's
    correlation is better than specified range (better means >= threshold
    or <= -threshold)
    """
    
    row += -2    # matches with what you see in the spreadsheet
    threshold = abs(threshold) / 100
    pear_corr = abs(df["pear_corr"][row])
    
    if pear_corr >= threshold:
        return True
    else:
        return False


# ================
def mir_stdev_test(
        df,
        row = 2,
        threshold = 50
        ):
# ================
    
    """This returns Bool after evaluating,
    in a DataFrame (likely produced by Reconciler), if the standard deviation
    in all miRNAs averages does not exceed threshold % of the averages.
    
    It expects these headers:
    mir_ave_b, mir_ave_p, mir_ave_l, mir_std_b,	mir_std_p, mir_std_l
    
    df = a pandas DataFrame
    row = the row where to run the test
    threshold = test will fail if the SD is beyond this average %
    
    If one of the averages is 0, False will be returned.
    """
    
    test = 0
    row += -2    # matches with what you see in the spreadsheet
    
    for g in ("b", "p", "l"):
        try:
            x = (
                (float(df["mir_std_" + g][row]) * 100)
                /  float(df["mir_ave_" + g][row])
                )
            if x >= threshold:
                test += 1
        except ZeroDivisionError:
            pass
            # if the average is *exactly* 0, we assume the SD is within range.
            # first, quantifications are all positives (so 0 can't be the
            # result of the average of non-zero elements), second the number
            # of allowed zeros in the averages is set elsewhere.
        except KeyError:
            print("Warning: key out of range. Use the row number as it appears"
                    " in the spreadsheed (e.g. first data row is 2)")
            raise
    
    if test > 0:
        return False
    else:
        return True


# ================
def lip_stdev_test(
        df = None,
        row = 0,
        threshold = 50
        ):
# ================
    
    """This returns Bool after evaluating,
    in a DataFrame (likely produced by Reconciler), if the standard deviation
    in all lipid averages does not exceed threshold % of the averages.
    
    It expects these headers:
    lip_ave_b, lip_ave_p, lip_ave_l, lip_std_b,	lip_std_p, lip_std_l
    
    df = a pandas DataFrame
    row = the row where to run the test
    threshold = test will fail if the SD is beyond this average %
    
    If one of the averages is 0, False will be returned.
    """
    
    test = 0
    row += -2    # matches with what you see in the spreadsheet
    
    for g in ("b", "p", "l"):
        try:
            x = (
                (float(df["lip_std_" + g][row]) * 100)
                /  float(df["lip_ave_" + g][row])
                )
            if x >= threshold:
                test += 1
        except ZeroDivisionError:
            pass
            # if the average is *exactly* 0, we assume the SD is within range.
            # first, quantifications are all positives (so 0 can't be the
            # result of the average of non-zero elements), second the number
            # of allowed zeros in the averages is set elsewhere.
        except KeyError:
            print("Warning: key out of range. Use the row number as it appears"
                    " in the spreadsheed (e.g. first data row is 2)")
            raise
    
    if test > 0:
        return False
    else:
        return True


# ============
def stdev_test(
        averages,
        standard_devs,
        threshold = 50
        ):
# ============
    
    """This returns Bool after evaluating,
    in supplied tuples/lists of values, if the standard deviation
    in all averages does not exceed threshold % of the averages.
    
    (!) It expects *matched* tuples/lists (the average in one must be
    associated with the standard deviations in the other one).
    *** Test is made row_index-wise. ***
    
    (!!) It also expects *same length* ! This is not checked but weird stuff
    is going to happen if that's not so.
    
    threshold = test will fail if the SD is beyond this average %
    
    If one of the averages is 0, False will be returned.
    """
    
    test = 0
    
    for v in range(len(averages)):
        try:
            x = (
                (float(standard_devs[v]) * 100)
                /  float(averages[v])
                )
            #print("*debug*stdev\t\tave\t\tx\t\tthreshold\n\t", standard_devs[v],
            #    "\t", averages[v], "\t", x, "\t", threshold)
            if x >= threshold:
                test += 1
        except ZeroDivisionError:
            pass
            # if the average is *exactly* 0, we assume the SD is within range.
            # first, quantifications are all positives (so 0 can't be the
            # result of the average of non-zero elements), second the number
            # of allowed zeros in the averages is set elsewhere.
    
    if test > 0:
        #print("*debug* TEST FAILED", test)
        return False
    else:
        #print("*debug* TEST PASSED", test)
        return True


# =====================================
def flush_to_file(string, outfile_name):
# =====================================
    
    """Checks whether outfile_name exists, then either writes or appends
    data from a string in memory."""
    
    if os.path.exists(outfile_name):
        with open(outfile_name, "a") as o:
            o.write(string)
    else:
        with open(outfile_name, "w") as o:
            o.write(string)
    


# =================================================
def list_permute_all(listlike1, listlike2):
# =================================================
    
    """This takes elements from two list-like objects, and returns two
    lists of equal length, each containing all possible index-to-index
    permutations.
    
    E.g.
    
    list1 = [10, 12, 14] # len = 3
    list2 = [1, 3, 5, 7] # len = 4 ; total permutations 3 x 4 = 12
    
    It returns:
    
    returned1 = [10, 10, 10, 10, 12, 12, 12, 12, 14, 14, 14, 14] # len = 12
    returned2 = [1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7] # len = 12
    """
    
    temp1, temp2 = [], []
    
    if len(listlike1) > len(listlike2):
        for elem in listlike2:
            for i in range(len(listlike1)):
                temp2.append(elem)
        for i in range(len(listlike2)):
            temp1 = temp1 + listlike1
    elif len(listlike1) < len(listlike2):
        for elem in listlike1:
            for i in range(len(listlike2)):
                temp1.append(elem)
        for i in range(len(listlike1)):
            temp2 = temp2 + listlike2
    elif len(listlike1) == len(listlike2):
        for i in range(len(listlike1)):
            temp1 = temp1 + listlike1
            temp2 = temp2 + listlike2
    
    return temp1, temp2


# =====================================================
def randomtest_pearson(dataset1, dataset2, points = 3):
# =====================================================
    
    """This function randomly picks <points> entries from
    dataset1 and dataset2, then performs pearson's correlation and
    returns (<pearson's value>, <p_value>)
    
    dataset1, dataset2: tuples or lists expected
    points: integer number of points being tested
    
    Tip:
    Since the pick is random with no particular distribution, dataset1
    and dataset2 should reflect the shape of the data being randomized.
    """
    
    if not isinstance(dataset1, (list, tuple)) or not isinstance(dataset2, (list, tuple)):
        raise TypeError("Datasets must be tuples or lists.")
    
    temp1 = [choice(dataset1) for i in range(points)]
    temp2 = [choice(dataset2) for i in range(points)]
    
    return pearsonr(temp1, temp2)

# =====================================================
def randomtest_spearman(dataset1, dataset2, points = 3):
# =====================================================
    
    """This function randomly picks <points> entries from
    dataset1 and dataset2, then performs pearson's correlation and
    returns (<pearson's value>, <p_value>)
    
    dataset1, dataset2: tuples or lists expected
    points: integer number of points being tested
    
    Tip:
    Since the pick is random with no particular distribution, dataset1
    and dataset2 should reflect the shape of the data being randomized.
    """
    
    if not isinstance(dataset1, (list, tuple)) or not isinstance(dataset2, (list, tuple)):
        raise TypeError("Datasets must be tuples or lists.")
    
    temp1 = [choice(dataset1) for i in range(points)]
    temp2 = [choice(dataset2) for i in range(points)]
    
    return spearmanr(temp1, temp2)


# =======================================================
def drop_rows_if_more_than_one_is_nan(dataframe, *header):
# =======================================================
    
    """This function checks at selected columns, row by row,
    and deletes the row from the dataframe if more than one value
    is NaN.
    Caution! If you used pd.DataFrame.replace() to change NaN to 0,
    the function will NOT treat 0 as NaN.
    
    Note: this does NOT modify input DataFrame.
    
    Returns: a modified copy of input DataFrame
    """
        
    dropme = [] # row indexes to drop at the end
    for r in range(df.shape[0]):
        hits = 0
        for h in list(header):
            if pd.isnull(df[h][r]):
                hits += 1
        if hits > 1:
            dropme.append(r)
        else:
            pass
    
    return dataframe.drop(dataframe.index[dropme])



# ==============================================
def count_repetitive_column_elements(df, header):
# ==============================================
    
    """Looks in a pandas dataframe <df> for all elements under <header>,
    and returns a dict holding the key : values for the ones that were found.
    
    Note: this dict can be then processed with pandas.DataFrame.from_dict to
    make a pandas df out of this dict.
    """
    
    mastermind = {}
    for id in df[header]:
        mastermind.setdefault(id, 0)
        mastermind[id] += 1

    return mastermind


# ======================================================
def write_dict_to_csv(dict, output_filename, sep = "\t"):
# ======================================================
    
    """This takes a dict and writes, as good as it can, a two-column csv
    out of the dict.
    Keys and values are passed to str().
    
    Beware of multidimensional dicts! This is primarily meant as a mean to
    write csv results produced by count_repetitive_column_elements() without
    having to do with pandas.DataFrame.from_dict, that usually complains
    about everything.
    
    You can specify your preferred separator with sep. Defaults to "\n"
    """
    
    holder = ""
    for key, val in dict.items():
        holder = holder + str(key) + sep + str(val) + "\n"
    with open(output_filename, "w") as o:
        o.write(holder)
        


# ==============================
def draw_scatterplot(
                array_1,    # x of dataset 1
                array_2,    # y of dataset 1
                array_3,    # x of dataset 2
                array_4,    # y of dataset 2
                array_5,    # x of dataset 3
                array_6,    # y of dataset 3
                polynomial_fit_grade = 1,   # 1 = line, > 1 = curve
                array_12_name = "first set",
                array_34_name = "second set",
                array_56_name = "third set",
                output_filename = "scatterplot.png"
                ):
# ==============================
    
    """This function aims at drawing a scatterplot of *three* different
    datasets, and also display a polynomial fit for it.
    It is designed to work with array-like data structures likewise arranged:
    
    dataset 1: x = l1, y = l2
    dataset 2: x = l3, y = l4
    dataset 3: x = l5, y = l6
    
    dataset names can be defined at function call (array_nn_name).
    """
        
    # changing datatypes from numpy to python ones for compatibility
    # float64s make Plotly go crazy
    l1 = [float(x) for x in array_1]
    l2 = [float(x) for x in array_2]
    l3 = [float(x) for x in array_3]
    l4 = [float(x) for x in array_4]
    l5 = [float(x) for x in array_5]
    l6 = [float(x) for x in array_6]
    
    trace1 = go.Scatter(
            x = l1,
            y = l2,
            mode = "markers",
            name = array_12_name
            )
    
    trace2 = go.Scatter(
            x = l3,
            y = l4,
            mode = "markers",
            name = array_34_name
            )
    
    trace3 = go.Scatter(
            x = l5,
            y = l6,
            mode = "markers",
            name = array_56_name
            )
    
    # -- fit calculation -- #
    
    lx = l1 + l3 + l5
    ly = l2 + l4 + l6
    
    # polyfit obtains the coefficients of a polynomial that fits the data
    # the order is given by polynomial_fit_grade
    z = polyfit(lx, ly, polynomial_fit_grade)
    
    # poly1d knows how to construct a polynomial from given array of
    # coefficients
    f = poly1d(z)
    
    # linspace gets evenly n spaced numbers between min and max
    if polynomial_fit_grade < 2:
        # we just get two points to draw a line for the fit
        x_new = linspace(min(lx), max(lx), 2)
    else:
        # we add more points to the fit if it is a curve
        x_new = linspace(min(lx), max(lx), len(lx))
        
    y_new = f(x_new)
    
    trace4 = go.Scatter(
            x = x_new, 
            y = y_new,
            mode = "lines",
            name = "Fit",
            )
    
    # // end of fit calculation // #
                
    data = [trace1, trace2, trace3, trace4]
    
    py.image.save_as({'data': data}, output_filename)
    
    # end of draw_scatterplot()


# modified version of draw_scatterplot to draw six datasets
# ==============================
def draw_6scatterplot(
                array_1,    # x of dataset 1
                array_2,    # y of dataset 1
                array_3,    # x of dataset 2
                array_4,    # y of dataset 2
                array_5,    # x of dataset 3
                array_6,    # y of dataset 3
                array_7,    # x of dataset 4
                array_8,    # y of dataset 4
                array_9,    # x of dataset 5
                array_a,    # y of dataset 5
                array_b,    # x of dataset 6
                array_c,    # y of dataset 6
                polynomial_fit_grade = 1,   # 1 = line, > 1 = curve
                array_12_name = "first set",
                array_34_name = "second set",
                array_56_name = "third set",
                array_78_name = "fourth set",
                array_9a_name = "fifth set",
                array_bc_name = "sixth set",          
                output_filename = "scatterplot.png"
                ):
# ==============================
    
    """This function aims at drawing a scatterplot of *six* different
    datasets, and also display a polynomial fit for it.
    It is designed to work with array-like data structures likewise arranged:
    
    dataset 1: x = l1, y = l2
    dataset 2: x = l3, y = l4
    dataset 3: x = l5, y = l6
    dataset 4: x = l7, y = l8
    dataset 5: x = l9, y = la
    dataset 6: x = lb, y = lc
    
    dataset names can be defined at function call (array_nn_name).
    """
        
    # changing datatypes from numpy to python ones for compatibility
    # float64s make Plotly go crazy
    l1 = [float(x) for x in array_1]
    l2 = [float(x) for x in array_2]
    l3 = [float(x) for x in array_3]
    l4 = [float(x) for x in array_4]
    l5 = [float(x) for x in array_5]
    l6 = [float(x) for x in array_6]
    l7 = [float(x) for x in array_7]
    l8 = [float(x) for x in array_8]
    l9 = [float(x) for x in array_9]
    la = [float(x) for x in array_a]
    lb = [float(x) for x in array_b]
    lc = [float(x) for x in array_c]
    
    trace1 = go.Scatter(
            x = l1,
            y = l2,
            mode = "markers",
            name = array_12_name
            )
    
    trace2 = go.Scatter(
            x = l3,
            y = l4,
            mode = "markers",
            name = array_34_name
            )
    
    trace3 = go.Scatter(
            x = l5,
            y = l6,
            mode = "markers",
            name = array_56_name
            )
    
    trace4 = go.Scatter(
            x = l7,
            y = l8,
            mode = "markers",
            name = array_78_name
            )
    
    trace5 = go.Scatter(
            x = l9,
            y = la,
            mode = "markers",
            name = array_9a_name
            )
    
    trace6 = go.Scatter(
            x = lb,
            y = lc,
            mode = "markers",
            name = array_bc_name
            )
    
    # -- fit calculation -- #
    
    lx = l1 + l3 + l5 + l7 + l9 + lb
    ly = l2 + l4 + l6 + l8 + la + lc
    
    # polyfit obtains the coefficients of a polynomial that fits the data
    # the order is given by polynomial_fit_grade
    z = polyfit(lx, ly, polynomial_fit_grade)
    
    # poly1d knows how to construct a polynomial from given array of
    # coefficients
    f = poly1d(z)
    
    # linspace gets evenly n spaced numbers between min and max
    if polynomial_fit_grade < 2:
        # we just get two points to draw a line for the fit
        x_new = linspace(min(lx), max(lx), 2)
    else:
        # we add more points to the fit if it is a curve
        x_new = linspace(min(lx), max(lx), len(lx))
        
    y_new = f(x_new)
    
    trace7 = go.Scatter(
            x = x_new, 
            y = y_new,
            mode = "lines",
            name = "Fit",
            )
    
    # // end of fit calculation // #
                
    data = [trace1, trace2, trace3, trace4, trace5, trace6, trace7]
    
    py.image.save_as({'data': data}, output_filename)
    
    # end of draw_6scatterplot()



# ==================================
def draw_chart_pear(array_1,
                    array_2,
                    dict,
                    pear_obj,
                    algorithm = "unk_algorithm",
                    sample = "unk_sample",
                    threshold = 0.0001
                    ):
# ==================================
    
    """Uses plot.ly and numpy/scipy to draw the regression fit to found values.
    
    It DOES expect to work on, in addition to two data arrays/lists, ALSO a
    extract_values() dict from which to pick some data to append the filename
    to.
    """
    
    if array_1 == None:
        raise ValueError("Empty array 1 to draw.")
    elif array_2 == None:
        raise ValueError("Empty array 2 to draw.")
    elif dict == None:
        raise ValueError("No extract_values() dict specified")
    elif pear_obj == None:
        raise ValueError("No Pearson's correlation object found "
                         "to read pval info from.")
    
    if pear_obj[1] <= threshold:
        
        trace1 = go.Scatter(
            x = array_1,
            y = array_2,
            mode = "markers"
            )
    
        z = polyfit(array_1, array_2, 1)
        f = poly1d(z)
        
        # linspace not needed for a line probably
        x_new = linspace(max(array_1), min(array_1), 2)
        y_new = f(x_new)
    
        trace2 = go.Scatter(
            x = x_new, 
            y = y_new,
            mode = "lines",
            name = "Fit",
            )
    
        data = [trace1, trace2]
        
        print("**debug** lipid prima", dict["curr_lipid_ID"])
        print("**debug** lipid dopo",
                dict["curr_lipid_ID"].replace("/", "_").replace(":", "¦")
                )
        
        filename = (str(pear_obj[1]) + "_" + sample
                    + "_" + algorithm + "_" + dict["curr_mirna_ID"]
                    + "_"
                    + dict["curr_lipid_ID"].replace("/", "_").replace(":", "¦")
                    + "_" + dict["curr_organ"]
                    + "_" + dict["curr_diet"]
                    + ".png")
    
        py.image.save_as({'data': data}, filename)
    else:
        pass



# ==============================
def draw_chart_permutations_pear(
                array_1,
                array_2,
                array_3,
                array_4,
                array_5,
                array_6,
                #dict,   # not supp curr
                pval,
                pcorr,
                mirna,
                organ,
                algorithm,
                lipid,
                sample,
                input_dataframe_name = "some_dataframe.csv",
                line = 0,
                array_12_name = "first set",
                array_34_name = "second set",
                array_56_name = "third set"
                ):
# ==============================
    
    """Uses plot.ly and numpy/scipy to draw the regression fit to found values.
    
    This is a modified version of draw_chart_pear() that is aimed at working
    with data coming from cooked mir + lip datasets (in the form of already
    found correlations), and at displaying the data in three different colors.
    """
    
    # don't let Nones as such. It will break.
    
    # there are lipid names out there that break the filesystem.
    lipid = lipid.replace(":", ".").replace("/", "-")
    
    # float64s make Plotly go crazy
    l1 = [float(x) for x in array_1]
    l2 = [float(x) for x in array_2]
    l3 = [float(x) for x in array_3]
    l4 = [float(x) for x in array_4]
    l5 = [float(x) for x in array_5]
    l6 = [float(x) for x in array_6]
    
    trace1 = go.Scatter(
            x = l1,
            y = l2,
            mode = "markers",
            name = array_12_name
            )
    
    trace2 = go.Scatter(
            x = l3,
            y = l4,
            mode = "markers",
            name = array_34_name
            )
    
    trace3 = go.Scatter(
            x = l5,
            y = l6,
            mode = "markers",
            name = array_56_name
            )
    
    lx = l1 + l3 + l5
    ly = l2 + l4 + l6
    
    z = polyfit(lx, ly, 1)
    f = poly1d(z)
    
    x_new = linspace(max(lx), min(ly), 2)
    y_new = f(x_new)
    
    trace4 = go.Scatter(
            x = x_new, 
            y = y_new,
            mode = "lines",
            name = "Fit",
            )
                
    data = [trace1, trace2, trace3, trace4]
    
    output_filename = (
            "row_" + str(line)
            + "_pval_" + str(pval)
            + "_pcorr_" + str(round(float(pcorr), 2))
            + "_" + mirna
            + "_" + algorithm
            + "_" + organ
            + "_lip_" + lipid
            + "_" + sample
            + ".png"
            )
    
    py.image.save_as({'data': data}, output_filename)
    
    print("\nWritten output image:\n" + output_filename + "\n")




# FIX! Make it general purpose
# =============================
def write_pear(outfile, string):
# =============================
    with open(outfile, "w") as o:
        o.write(string)


# ========== DEBUG FUNCTION ======================
def debug_check_values_lenght(dictionary, logfile):
# ================================================
    
    """This is a debug function that checks if the number of expression
    levels pulled are of the desired length.
    """

    if (len(dictionary["mirna_expression_levels_b_desV2"]) == 3 and
        len(dictionary["mirna_expression_levels_p_desV2"]) == 3 and
        len(dictionary["mirna_expression_levels_l_desV2"]) == 3 and
        len(dictionary["mirna_expression_levels_b_UQ"]) == 3 and
        len(dictionary["mirna_expression_levels_p_UQ"]) == 3 and
        len(dictionary["mirna_expression_levels_l_UQ"]) == 3):
        pass
    else:
        print("**debug** MIRNA TEST FAILED",
                dictionary["curr_diet"],
                dictionary["curr_organ"],
                dictionary["curr_mirna_ID"],
                dictionary["curr_lipid_ID"]
                )

        logfile = (logfile
                + get_current_time()
                + "**debug** MIRNA TEST FAILED "
                + dictionary["curr_diet"]
                + " ," + dictionary["curr_organ"]
                + " ," + dictionary["curr_mirna_ID"]
                + " ," + dictionary["curr_lipid_ID"] + "\n"
                + "mirna_expression_levels_b_desV2 length: "
                + str(len(dictionary["mirna_expression_levels_b_desV2"]))
                + str(dictionary["mirna_expression_levels_b_desV2"]) + "\n"
                + "mirna_expression_levels_p_desV2 length: "
                + str(len(dictionary["mirna_expression_levels_p_desV2"]))
                + str(dictionary["mirna_expression_levels_p_desV2"]) + "\n"
                + "mirna_expression_levels_l_desV2 length: "
                + str(len(dictionary["mirna_expression_levels_l_desV2"]))
                + str(dictionary["mirna_expression_levels_l_desV2"]) + "\n"
                + "mirna_expression_levels_b_UQ length: "
                + str(len(dictionary["mirna_expression_levels_b_UQ"]))
                + str(dictionary["mirna_expression_levels_b_UQ"]) + "\n"
                + "mirna_expression_levels_p_UQ length: "
                + str(len(dictionary["mirna_expression_levels_p_UQ"]))
                + str(dictionary["mirna_expression_levels_p_UQ"]) + "\n"
                + "mirna_expression_levels_l_UQ length: "
                + str(len(dictionary["mirna_expression_levels_l_UQ"]))
                + str(dictionary["mirna_expression_levels_l_UQ"]) + "\n***\n"
                )
        
        return logfile


    if (len(dictionary["lipid_amount_b_pla"]) == 6 and
        len(dictionary["lipid_amount_p_pla"]) == 6 and
        len(dictionary["lipid_amount_l_pla"]) == 6 and
        len(dictionary["lipid_amount_b_aor"]) == 6 and
        len(dictionary["lipid_amount_p_aor"]) == 6 and
        len(dictionary["lipid_amount_l_aor"]) == 6 and
        len(dictionary["lipid_amount_b_liv"]) == 6 and
        len(dictionary["lipid_amount_p_liv"]) == 6 and
        len(dictionary["lipid_amount_l_liv"]) == 6):
        #print("**debug** lipi lengths pass test", dictionary["curr_diet"], dictionary["curr_organ"])
        pass
    else:
        print("**debug** LIPI TEST FAILED", 
                dictionary["curr_diet"],
                dictionary["curr_organ"],
                dictionary["curr_mirna_ID"],
                dictionary["curr_lipid_ID"]
                )
        logfile = (logfile
                + get_current_time()
                + "**debug** LIPI TEST FAILED "
                + dictionary["curr_diet"]
                + " ," + dictionary["curr_organ"]
                + " ," + dictionary["curr_mirna_ID"]
                + " ," + dictionary["curr_lipid_ID"] + "\n"
                + "lipid_amount_b_pla length, content: "
                + str(len(dictionary["lipid_amount_b_pla"]))
                + str(dictionary["lipid_amount_b_pla"]) + "\n"
                + "lipid_amount_p_pla llength, content: "
                + str(len(dictionary["lipid_amount_p_pla"]))
                + str(dictionary["lipid_amount_p_pla"]) + "\n"
                + "lipid_amount_l_pla length, content: "
                + str(len(dictionary["lipid_amount_l_pla"]))
                + str(dictionary["lipid_amount_l_pla"]) + "\n"
                + "lipid_amount_b_aor length, content: "
                + str(len(dictionary["lipid_amount_b_aor"]))
                + str(dictionary["lipid_amount_b_aor"]) + "\n"
                + "lipid_amount_p_aor length, content: "
                + str(len(dictionary["lipid_amount_p_aor"]))
                + str(dictionary["lipid_amount_p_aor"]) + "\n"
                + "lipid_amount_l_aor length, content: "
                + str(len(dictionary["lipid_amount_l_aor"]))
                + str(dictionary["lipid_amount_l_aor"]) + "\n"      
                + "lipid_amount_b_liv length, content: "
                + str(len(dictionary["lipid_amount_b_liv"]))
                + str(dictionary["lipid_amount_b_liv"]) + "\n"
                + "lipid_amount_p_liv length, content: "
                + str(len(dictionary["lipid_amount_p_liv"]))
                + str(dictionary["lipid_amount_p_liv"]) + "\n"
                + "lipid_amount_l_liv length, content: "
                + str(len(dictionary["lipid_amount_l_liv"]))
                + str(dictionary["lipid_amount_l_liv"]) + "\n***\n"
                )
        
        return logfile


# =========================    
def get_lipid_class(string):
# =========================
    
    try:
        match = re.search(r"(CE|Cer|DAG|FC|Gb3|Glc|Lac|LPC|LPE|LPI|PA|PC \d|PC O|PC P|PE \d|PE P|PE O|PG|PI|PS|SM|TAG)", string)
        return match.group(1)
    except:
        raise TypeError("Could not determine lipid class for:", string)


# ==========================
def get_lipid_class2(string):
# ==========================
    
    """
    Improved version of former function, now it groups PC and PE
    lipid classes regardless of the digits that make up lipid name.
    
    Returns: <string> of lipid class
    """
    
    # when working with whole classes, there's no PC or PE followed by digits.
    if string in ("PC", "PE"):
        return string
    
    try:
        match = re.search(r"(CE|Cer|DAG|FC|Gb3|Glc|Lac|LPC|LPE|LPI|PA|PC \d|PC O|PC P|PE \d|PE P|PE O|PG|PI|PS|SM|TAG)", string)
        matched = match.group(1)
        
        # grouping PC \d and PE \d classes
        if re.search(r"PC \d", matched):
            return "PC"
        elif re.search(r"PE \d", matched):
            return "PE"
        else:
            return matched
        
    except:
        fatty_acids = [
            "17:1", "13:2", "14:1", "15:0", "20:3", "13:1",
            "20:1", "19:0", "17:0", "22:2", "19:1", "14:0",
            "22:1", "21:1", "16:1", "18:1", "18:2", "20:4",
            "18:0", "22:3", "22:4", "23:0", "16:0", "22:6",
            "19:2", "23:1", "25:1", "18:3", "12:0", "22:5",
            "20:5", "20:2", "20:0", "24:1", "15:1", "23:2",
            "24:0", "26:0", "15:2", "26:1", "21:0", "22:0",
            "24:5", "25:2", "24:2"
        ]
        if string in fatty_acids:
            return string
        else:
            raise TypeError("Could not determine lipid class for:", string)


# ===================
match_lipi_ids_func = re.compile(r"(..)_(.*)_(.*?)_([^_]*)")
def lipi_ids(string):
# ===================
    
    """It reads and parses the headers of adjusted lipidomics table.
    These look like: AO05_LDLR_1f_chow
    It returns a tuple containing translated and ordered elements, likewise:
    (sample, genotype, diet, animal)
    """
    
    betternames = {
        "LI" : "liv",
        "AO" : "aor",
        "PL" : "pla",
        "PCSK9" : "p",
        "BL6" : "b",
        "LDLR" : "l",
        "WD" : "wd",
        "chow" : "cd"
    }
    
    try:
        match = match_lipi_ids_func.search(string)
        sample = betternames[match.group(1)] # returned 0
        genotype = betternames[match.group(2)]   # returned 1
        animal = match.group(3) # returned 3
        diet = betternames[match.group(4)]  # returned 2
        return (sample, genotype, diet, animal)
    except:
        raise TypeError("Something like AO_LDLR_1f_chow expected, instead:", string)


# ====================
match_mirnaids_func = re.compile(r"(.)(.[^_]*)_(.[^_]*)_(.[^_]*)_(.[^_]*)")
def mirna_ids(string):
# ====================

    """It reads and parses the headers of adjusted mirnomics table.
    These look like: l6_ile_wd_desV2
    It returns a tuple containing translated and ordered elements, likewise:
    (organ, genotype, diet, animal, algorithm)
    """
    
    betternames = {
        "chow" : "cd",
        "wd" : "wd"
    }
    
    try:
        match = match_mirnaids_func.search(string)
        genotype = match.group(1)   # returned 1
        animal = match.group(2) # returned 4
        organ = match.group(3)  #returned 0 | liv aor bra wat jej duo ile
        diet = betternames[match.group(4)]  # returned 2
        algorithm = match.group(5)  #returned 5 | either UQ or desV2
        return (organ, genotype, diet, animal, algorithm)
    except:
        raise TypeError("Something like b3m_aor_chow_desV2 expected.")


# =============
def gen_cond():
# =============
    
    """Conditions (mirna) generator.
    Note: genotype is not generated!
    """
    
    # our playground
    mirna_organs = ("aor", "liv", "bra", "duo", "jej", "ile", "wat")
    #genotypes = ("b", "p", "l")
    diets = ("cd", "wd")
    #algorithms = ("desV2", "UQ")
    
    for d in diets:
        for o in mirna_organs:
            yield (o, d)    # NOT the same pattern as mirna_ids()



# =======================
def gen_fake_lipid_cond():
# ======================
    
    """Generates fake headers for lipidomics datatable compatible with
    lipi_ids()
    """
    
    # our playground
    samples = ("LI", "AO", "PL")
    genotypes = ("PCSK9", "BL6", "LDLR")
    diets = ("chow", "WD")
    
    for s in samples:
        for g in genotypes:
            for d in diets:
                yield (s + "_" + g + "_00_" + d)


# ============================
def gen_fake_lipid_cond_full(
        biological_replicates = 6
        ):
# ============================
    
    """Generates fake headers for lipidomics datatable compatible with
    lipi_ids()
    
    It also generates desired number of biological replicates.
    """
    
    # our playground
    samples = ("LI", "AO", "PL")
    genotypes = ("PCSK9", "BL6", "LDLR")
    diets = ("chow", "WD")
    
    for s in samples:
        for g in genotypes:
            for d in diets:
                for i in range(biological_replicates):
                    yield (s + "_" + g + "_" + str(i) + "_" + d)


# ==============================
def gen_fake_lipid_cond_full_sex(
        biological_replicates = 6
        ):
# ==============================
    
    """Generates fake headers for lipidomics datatable compatible with
    lipi_ids()
    
    It also generates desired number of biological replicates.
    """
    
    # our playground
    samples = ("LI", "AO", "PL")
    genotypes = ("PCSK9", "BL6", "LDLR")
    diets = ("chow", "WD")
    
    males = int(biological_replicates / 2)
    
    for s in samples:
        for g in genotypes:
            for d in diets:
                for i in range(biological_replicates):
                    if i < males:
                        yield (s + "_" + g + "_" + str(i) + "m_" + d)
                    else:
                        yield (s + "_" + g + "_" + str(i) + "f_" + d)


# =======================
def gen_fake_mirna_cond():
# ======================
    
    """Generates fake headers for mirnomics datatable compatible with
    mirna_ids()
    """
    
    # our playground
    samples = ("aor", "bra", "wat", "duo", "jej", "ile", "liv")
    genotypes = ("b", "l", "p")
    diets = ("chow", "wd")
    algorithms = ("desV2", "UQ")
    
    for a in algorithms:
        for d in diets:
            for s in samples:
                for g in genotypes:
                    yield (g + "00_" + s + "_" + d + "_" + a)


# ===============================
def gen_fake_mirna_cond_full(
        biological_replicates = 3
        ):
# ===============================
    
    """Generates fake headers for mirnomics datatable compatible with
    mirna_ids()
    
    It also generates desired number of biological replicates.
    """
    
    # our playground
    samples = ("aor", "bra", "wat", "duo", "jej", "ile", "liv")
    genotypes = ("b", "p", "l")
    diets = ("chow", "wd")
    algorithms = ("desV2", "UQ")
    
    for a in algorithms:
        for d in diets:
            for s in samples:
                for g in genotypes:
                    for i in range(biological_replicates):
                       yield (g + str(i) + "_" + s + "_" + d + "_" + a)


# ==============================
def gen_fake_mirna_cond_full_sex(
        biological_replicates = 3
        ):
# ==============================
    
    """Generates fake headers for mirnomics datatable compatible with
    mirna_ids(); the number of animals is pushed towards higher values
    so that it does not overlap with fake lipid ones.
    Note: no sex is evaluated as male by Reconciler --statsfriendly.
        
    It also generates desired number of biological replicates.
    """
    
    # our playground
    samples = ("aor", "bra", "wat", "duo", "jej", "ile", "liv")
    genotypes = ("b", "p", "l")
    diets = ("chow", "wd")
    algorithms = ("desV2", "UQ")
    
    for a in algorithms:
        for d in diets:
            for s in samples:
                for g in genotypes:
                    for i in range(biological_replicates):
                       yield (g + str(i + 5) + "_" + s + "_" + d + "_" + a)




# =====================
def get_current_time():
# =====================
    
    """ This returns a properly formatted string with local time.
    """
    
    now = time.localtime()
    
    weekdays = {
        "0" : "Monday",
        "1" : "Tuesday",
        "2" : "Wednesday",
        "3" : "Thursday",
        "4" : "Friday",
        "5" : "Saturday",
        "6" : "Sunday"
    }
    
    right_now = (weekdays[str(now[6])]
            + ", " + str(now[2])
            + "/" + str(now[1])
            + "/" + str(now[0])
            + " - " + str(now[3])
            + ":" + str(now[4])
            + ":" + str(now[5])
            )
    
    return right_now


# ========================
def strip_mirname(string):
# ========================
    
    """Takes:
    mmu-miR-212-5p
    out of:
    mmu-miR-212-5p_ID=MIM0017053_MI0000696_mmu-mir-212
    """
    
    # 07/10/24 hack to work with anon data
    if string.startswith("xxx"):
        return string[4:]
        
    try:
        match = re.search(r"(.[^_]*)_.*", string)
        return match.group(1)
    except:
        raise TypeError(
            "Something like "
            "mmu-miR-212-5p_ID=MIM0017053_MI0000696_mmu-mir-212 expected. "
            f"{string} got instead."
            )
    

# FIXFIXFIX: I have the feeling this can be done better.
# save cycles pulling values from df by giving *direct* header names
# instead of cycling all of them and doing the job for just the right ones
# that are the VAST minority of the whole.



# ====================
def add_averages(dict):
# ====================
    
    """ This function works on dicts generated by extract_values().
    This function calculates the average for mirnas expression levels as
    well as lipid amounts and stores them in the dictionary, which is returned.
    
    Please note that * THE INPUT IS RETURNED MODIFIED * and this is the
    intended behavior. Next time I'll consider doing a class with proper
    methods ._.'
    """
    
    dict["b_desV2_ave"] = average(dict["mirna_expression_levels_b_desV2"])
    dict["p_desV2_ave"] = average(dict["mirna_expression_levels_p_desV2"])
    dict["l_desV2_ave"] = average(dict["mirna_expression_levels_l_desV2"])
    dict["b_UQ_ave"] = average(dict["mirna_expression_levels_b_UQ"])
    dict["p_UQ_ave"] = average(dict["mirna_expression_levels_p_UQ"])
    dict["l_UQ_ave"] = average(dict["mirna_expression_levels_l_UQ"])
    
    dict["b_lip_aor_ave"] = average(dict["lipid_amount_b_aor"])
    dict["p_lip_aor_ave"] = average(dict["lipid_amount_p_aor"])
    dict["l_lip_aor_ave"] = average(dict["lipid_amount_l_aor"])
    dict["b_lip_liv_ave"] = average(dict["lipid_amount_b_liv"])
    dict["p_lip_liv_ave"] = average(dict["lipid_amount_p_liv"])
    dict["l_lip_liv_ave"] = average(dict["lipid_amount_l_liv"])
    dict["b_lip_pla_ave"] = average(dict["lipid_amount_b_pla"])
    dict["p_lip_pla_ave"] = average(dict["lipid_amount_p_pla"])
    dict["l_lip_pla_ave"] = average(dict["lipid_amount_l_pla"])
    
    # this does not return anything. The object is itself modified.





# =============================================================================
#                           * end of functions *
# =============================================================================