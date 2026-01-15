import matplotlib.pyplot as plt
import numpy as np

from ddg_func_graph import ddg_scores, func_scores 


def create_csv(inp, inp2):
    """Creates graph of functionality and ddg

    Args:
        inp (file path): location of the ddg data
        inp2 (file path): location of the functionality data
    """
    
    # read in the ddg scores csv
    aa_res_list, roll_aa_ddg_list, max_res = ddg_scores(inp, False)
    
    # read up functionality scores csv
    aa_res_list_f, roll_aa_func_list = func_scores(inp2, False)

    # Sort to only include functionality of high ddg
    aa_res_list, true_func_list = sorted_func(aa_res_list, roll_aa_ddg_list, aa_res_list_f, roll_aa_func_list)

    # create 2 lists of func. Above and below 0.5:
    aa_res_list_large, func_list_large, aa_res_list_small, func_list_small = func_split(aa_res_list, true_func_list)


    # font sizes
    titleSize = 20 
    labelSize = 14
    descSize = 10

    # create the graph
    fig, ax = plt.subplots(figsize=[8.4, 4.8])
    marker_dict = {"none":"none","point":".","circle":"o","triangle down":"v","octagon":"8","square":"s","star":"*"}
    ax.set_title("Functionality of Stable Residues", fontsize=titleSize)
    ax.set_xlabel("Residue #", fontsize=labelSize)


    #create the functionality graph
    ax.set_ylabel("Fuctionality Scores", fontsize=labelSize)
    marker_ = marker_dict["point"]
    ax.scatter(aa_res_list_large, y=func_list_large, c="red", marker=marker_, s=36.0, label="High Functionality")
    ax.scatter(aa_res_list_small, y=func_list_small, c="blue", marker=marker_, s=36.0, label="Low Functionality")

    # annotations
    ax.axvspan(1, 35, color='#FFA3A3', alpha=0.3)
    ax.annotate('Ser\nRich',
            xy=(2, -0.1),
            horizontalalignment='left', verticalalignment='top',
            fontsize=descSize)

    ax.axvspan(48, 296, color='#7CAED1', alpha=0.3)
    ax.axvspan(296, 338, color='#7CAED1', alpha=0.3)
    ax.annotate('Catalytic\nDomain',xy=(120, -0.1),horizontalalignment='left', verticalalignment='top',fontsize=descSize)
    ax.annotate('Hinge\nRegion',xy=(298, -0.1),horizontalalignment='left', verticalalignment='top',fontsize=8)

    
    ax.axvspan(338, 381, color='#C770CE', alpha=0.3)
    ax.axvspan(381, 394, color='#C770CE', alpha=0.3)
    ax.axvspan(394, 412, color='#C770CE', alpha=0.3)
    ax.annotate('Linker\nRegion',xy=(357, -0.1),horizontalalignment='left', verticalalignment='top',fontsize=descSize)
    ax.annotate('L1',xy=(352, 0.0),horizontalalignment='left', verticalalignment='top',fontsize=8)
    ax.annotate('L2',xy=(382, 0.0),horizontalalignment='left', verticalalignment='top',fontsize=8)
    ax.annotate('L3',xy=(397, 0.0),horizontalalignment='left', verticalalignment='top',fontsize=8)

    ax.axvspan(412, 644, color='#B79900', alpha=0.3)
    ax.annotate('Regulatory\nDomain',xy=(500, -0.1),horizontalalignment='left', verticalalignment='top',fontsize=descSize)
    
    ax.plot([-10, 337], [0.25,0.25], marker = 'o', color="#2c40b0",markersize=0)
    ax.plot([338, max_res+10], [0.70,0.70], marker = 'o', color="#2c40b0",markersize=0)


    
    #finilize and show graph
    ax.legend(fontsize=labelSize, loc="center right")
    ax.set_xbound(-5, max_res+5)
    ax.set_ybound(-0.2,1)
    
    plt.style.context("seaborn-talk")
    fig.savefig("graphs/highDDG_lowFunc.png", bbox_inches='tight')
    plt.show()


def sorted_func(aa_res_list, roll_aa_ddg_list, aa_res_list_f, roll_aa_func_list):

    # go through each ddg, if it is below -0.1, remove it from aa_func lists
    for ind, ddg in enumerate(roll_aa_ddg_list):
        if ddg < -0.1:

            # get res number
            res_num = aa_res_list[ind]

            # get index number in func list
            func_ind = aa_res_list_f.index(res_num)

            # remove that index from
            del aa_res_list_f[func_ind]
            del roll_aa_func_list[func_ind]

    return aa_res_list_f, roll_aa_func_list


def func_split(aa_res_list, true_func_list):
    aa_res_list_large = []
    func_list_large = []
    aa_res_list_small = []
    func_list_small = [] 

    for ind in range(len(aa_res_list)):
        if(aa_res_list[ind] >= 338):
            start = ind
            break
        if true_func_list[ind] > 0.25:
            aa_res_list_large.append(aa_res_list[ind])
            func_list_large.append(true_func_list[ind])
        else:
            aa_res_list_small.append(aa_res_list[ind])
            func_list_small.append(true_func_list[ind])
        
    for ind2 in range(len(aa_res_list[start:])):
        ind = ind2 + start
        if true_func_list[ind] > 0.70:
            aa_res_list_large.append(aa_res_list[ind])
            func_list_large.append(true_func_list[ind])
        else:
            aa_res_list_small.append(aa_res_list[ind])
            func_list_small.append(true_func_list[ind])

    return aa_res_list_large, func_list_large, aa_res_list_small, func_list_small



if __name__ == "__main__":
    inp = r"D:\Research Lab\ddGun - Functionality Analysis\ddgun output\ddg_op.txt"
    inp2 = r"D:\Research Lab\ddGun - Functionality Analysis\ddgun output\f25WT_t1_simple_aa_floored.csv"
    
    create_csv(inp, inp2)