import matplotlib.pyplot as plt
import numpy as np

from ddg_func_graph import create_ddg_func_graph 
from only_func_graph import create_func_graph 
from only_ddg_graph import create_ddg_graph 



if __name__ == "__main__":
    inp = r"D:\Research Lab\ddGun - Functionality Analysis\ddgun output\ddg_op.txt"                       #F for PC, D for lap
    inp2 = r"D:\Research Lab\ddGun - Functionality Analysis\ddgun output\f25WT_t1_simple_aa_floored.csv"  #F for PC, D for lap
    typ = "func"  # can be "all", "func", "ddg"
    roll = True  # can be True or False
    rg = ()      # can be empty (then entire), or 2 floats


    
    if(typ == "all"):
        create_ddg_func_graph(inp, inp2, roll, rg)
    elif(typ == "func"):
        create_func_graph(inp, inp2, roll, rg)
    elif(typ == "ddg"):
        create_ddg_graph(inp, inp2, roll, x_range = ())