import csv

def create_init(pangolearn_new_v,pango_version,outfile):

    with open(outfile,"w") as fw:
        fw.write(f'''_program = "pangoLEARN"
__version__ = "{pangolearn_new_v}"
PANGO_VERSION = "{pango_version}"

__all__ = [
    "training",
    "utils"]

from pangoLEARN import *

''')


def get_dict(in_csv,name_column,data_column):
    this_dict = {}
    with open(in_csv,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            this_dict[row[name_column]] = row[data_column]
    return this_dict
