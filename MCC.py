import numpy as np
from re import finditer
from os import listdir
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
from scipy.stats import ttest_ind

def calc_MCC(prediction, actual, name, out_number, alg):
    TP = 0
    TN = 0
    FP = 0
    FN = 0
    with open("{}_{}_{}_colormap.col".format(name, out_number, alg), 'w+') as f:

        for i, (a, b) in enumerate(zip(prediction, actual)):
            #true
            if a == b:
                f.write("{}:green ".format(i + 1)) #forna is 1-indexed...
                #negative
                if b == -1:
                    TN += 1
                #positive
                else:
                    TP += 1

            #false
            else:
                f.write("{}:red ".format(i + 1))
                #positive
                if b == -1:
                    FP += 1
                #negative
                else:
                    FN += 1

    MCC = ((TP * TN) - (FP * FN)) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

    return MCC

#I represent structures as a pairing vector
def parse_dot_bracket(input):
    output = np.full(len(input), -1)
    #I'm not sure this is the most efficent way to do this, but I'm lazy.
    more = True
    while more:
        more = False

        #finds matched parenthesis
        for x in finditer(r"\([^()]*\)", input):
            more = True
            output[x.start()] = x.end()-1
            output[x.end()-1] = x.start()

            #its recursive...
            input=input[0:x.start()] + "." + input[x.start()+1:x.end()-1] + "." + input[x.end():]

    return output

def make_fasta_db(name, out_number, alg, seq, struct):
    with open("{}_{}_{}.fasta".format(name, out_number, alg), 'w+') as f:
        f.write(">"+name+"\n")
        f.write(seq+"\n")
        f.write(struct+"\n")

def main(rank):

    step = 4
    out_number = rank
    if rank == 0:
        rank = 1
        step = 1

    #read all files in the directory, extract dot-brackets and calculate MCC
    data_fold = []
    data_roll = []
    data_diff = []
    actual_energy = []
    RNAfold_energy = []
    rollout_energy = []
    runtime = []
    actual_MFE = []
    RNAfold_MFE = []
    rollout_MFE = []
    names = []
    seq_lens = []
    families = []
    for fname in listdir():
        if fname.split(".")[-1] == "csv":
            with open(fname, "r") as f:
                lines = f.readlines()
                for l in lines[rank::step]:
                    l = l.split(",")
                    if l[0] == 'NONE':
                        continue
                    name = l[0]
                    seq = l[1]
                    seq_lens.append(len(seq))
                    actual = l[2]
                    RNAfold = l[3]
                    rollout = l[4]
                    family = name.split("_")[0]
                    actual_energy.append(float(l[5]))
                    RNAfold_energy.append(float(l[6]))
                    rollout_energy.append(float(l[7]))
                    runtime.append(float(l[15]))
                    actual_MFE.append(float(l[25]))
                    RNAfold_MFE.append(float(l[26]))
                    rollout_MFE.append(float(l[27]))
                    
                    make_fasta_db(name, out_number, "actual", seq, actual)
                    make_fasta_db(name, out_number, "RNAfold", seq, RNAfold)
                    make_fasta_db(name, out_number, "rollout", seq, rollout)
                    actual = parse_dot_bracket(actual)            
                    rollout = parse_dot_bracket(rollout)
                    RNAfold = parse_dot_bracket(RNAfold)
                    acc_fold = calc_MCC(RNAfold, actual, name, out_number, "RNAfold")
                    acc_roll = calc_MCC(rollout, actual, name, out_number, "rollout")
                    data_fold.append(acc_fold)
                    data_roll.append(acc_roll)
                    data_diff.append(acc_roll - acc_fold)
                    names.append(name)
                    families.append(family)
                    #print(name, acc_fold, acc_roll, acc_roll-acc_fold)

    all_data = np.array([data_fold, data_roll, data_diff, actual_energy, RNAfold_energy, rollout_energy, seq_lens, runtime, actual_MFE, RNAfold_MFE, rollout_MFE], dtype=float)
    DATA_FOLD = 0
    DATA_ROLL = 1
    DATA_DIFF = 2
    ACTUAL_FOLDABILITY = 3    
    RNAFOLD_FOLDABILITY = 4
    ROLLOUT_FOLDABILITY = 5
    SEQ_LEN = 6
    RUNTIME = 7
    ACTUAL_MFE = 8
    RNAFOLD_MFE = 9
    ROLLOUT_MFE = 10
    sort = all_data[DATA_DIFF].argsort()
    all_data = all_data[:,sort]
    names = np.array(names, dtype=str)[sort]
    families = np.array(families)[sort]

    statistic, pvalue = ttest_ind(all_data[DATA_FOLD], all_data[DATA_ROLL])

    print("{},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f}".format(out_number, np.mean(all_data[ROLLOUT_FOLDABILITY]), np.mean(all_data[DATA_ROLL]), np.median(all_data[DATA_ROLL]), np.mean(all_data[DATA_DIFF]), np.std(all_data[DATA_DIFF]), pvalue))

    with open("scores_{}.dat".format(out_number), 'w') as f:
        for n, d in zip(names, all_data.T):
            f.write("{} {} {} {}\n".format(n, d[DATA_FOLD], d[DATA_ROLL], d[DATA_DIFF]))


    mpl.rcParams.update({'font.size': 22})
    mpl.rc('text', usetex=True)
    mpl.rc('xtick', labelsize=18)
    mpl.rc('ytick', labelsize=18)

    fnames = sorted(list(set(families)))
    best_RNAfold = all_data[DATA_FOLD] > 0.9
    best_rollout = all_data[DATA_ROLL] > 0.9

    #comparison scatterplot
    fig, ax = plt.subplots()
    for fname in fnames:
        ax.scatter(all_data[DATA_FOLD][families == fname], all_data[DATA_ROLL][families == fname], label=fname, alpha=0.4)
    ax.set_xlabel("RNAfold MCC")
    ax.set_ylabel("ExpertRNA MCC")
    ax.set_xlim(-1, 1.1)
    ax.legend(fontsize=14)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)    
    plt.tight_layout()
    plt.savefig("scatter_comparison_{}.pdf".format(out_number))
    plt.close()

    #histograms!
    #overlay
    fig, ax = plt.subplots()
    ax.hist(data_fold, bins=50, range=(-1, 1), alpha=0.4, label="RNAfold")
    ax.hist(data_roll, bins=50, range=(-1, 1), alpha=0.4, label="ExpertRNA")
    ax.set_xlabel("Matthews Correlation Coefficient")
    ax.set_ylabel("Number")
    ax.set_xlim(-1, 1.2)
    ax.legend(fontsize=14)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig("comparison_{}.pdf".format(out_number))
    plt.close()

    #difference
    fig, ax = plt.subplots()
    ax.hist(data_diff, range=(-1,1), bins=30)
    ax.set_xlabel("MCC difference (ExpertRNA-RNAfold)")
    ax.set_ylabel("Number")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(-1, 1.2)
    plt.tight_layout()
    plt.savefig("diff_{}.pdf".format(out_number))
    plt.close()

    #runtime
    fig, ax = plt.subplots()
    def func(x, a, b, c):
        return a * np.power(2, b*x) + c
    def func2(x, a, b, c):
        return a*np.power(x, 2) + b * x + c
    #popt, pcov = curve_fit(func, all_data[SEQ_LEN], all_data[RUNTIME])
    popt2, pcov2 = curve_fit(func2, all_data[SEQ_LEN], all_data[RUNTIME])
    mi = min(all_data[SEQ_LEN])
    ma = max(all_data[SEQ_LEN])
    modelX = np.linspace(mi, ma, 100)
    modelY = func2(modelX, *popt2)
    ax.scatter(all_data[SEQ_LEN], all_data[RUNTIME])
    ax.plot(modelX, modelY, label = r"quadratic fit ${:.2f}x^{{2}} + {:.2f}x + {:.2f}$".format(popt2[0], popt2[1], popt2[2]))
    ax.set_xlabel("Sequence Length")
    ax.set_ylabel("Runtime (s)")
    ax.legend(fontsize=14)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig("runtime_{}.pdf".format(out_number))
    plt.close()

    #actual vs rollout MFE
    fig, ax = plt.subplots()
    line = np.linspace(-0.6, 0, 20)
    for fname in fnames:
        ax.scatter(all_data[ACTUAL_MFE][families == fname] / all_data[SEQ_LEN][families == fname], all_data[ROLLOUT_MFE][families == fname] / all_data[SEQ_LEN][families == fname], label=fname, alpha=0.4)
    ax.scatter(all_data[ACTUAL_MFE][best_rollout] / all_data[SEQ_LEN][best_rollout], all_data[ROLLOUT_MFE][best_rollout] / all_data[SEQ_LEN][best_rollout], c='red', s=8, label="Best ExpertRNA")
    ax.scatter(all_data[ACTUAL_MFE][best_RNAfold] / all_data[SEQ_LEN][best_RNAfold], all_data[ROLLOUT_MFE][best_RNAfold] / all_data[SEQ_LEN][best_RNAfold], c='cyan', s=4, label="Best RNAfold")
    ax.plot(line, line, c='k', linewidth=0.75)
    ax.set_xlabel("Actual MFE (normalized)")
    ax.set_ylabel("ExpertRNA MFE (normalized)")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axis('scaled')
    leg = ax.legend(fontsize=14, bbox_to_anchor=(1.04, 1))
    plt.tight_layout()
    plt.savefig("actual_roll_MFE_comparison_{}.pdf".format(out_number), bbox_extra_artists=[leg], bbox_inches='tight')
    plt.close()

    #actual vs RNAfold MFE
    fig, ax = plt.subplots()
    for fname in fnames:
        ax.scatter(all_data[ACTUAL_MFE][families == fname] / all_data[SEQ_LEN][families == fname], all_data[RNAFOLD_MFE][families == fname] / all_data[SEQ_LEN][families == fname], label=fname, alpha=0.4)
    ax.scatter(all_data[ACTUAL_MFE][best_rollout] / all_data[SEQ_LEN][best_rollout], all_data[RNAFOLD_MFE][best_rollout] / all_data[SEQ_LEN][best_rollout], c='red', s=8, label="Best ExpertRNA")
    ax.scatter(all_data[ACTUAL_MFE][best_RNAfold] / all_data[SEQ_LEN][best_RNAfold], all_data[RNAFOLD_MFE][best_RNAfold] / all_data[SEQ_LEN][best_RNAfold], c='cyan', s=4, label="Best RNAfold")
    ax.plot(line, line, c='k', linewidth=0.75)
    ax.set_xlabel("Actual MFE (normalized)")
    ax.set_ylabel("RNAfold MFE (normalized)")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axis('scaled')
    leg = ax.legend(fontsize=14, bbox_to_anchor=(1.04, 1))
    plt.tight_layout()
    plt.savefig("actual_RNAfold_MFE_comparison_{}.pdf".format(out_number), bbox_extra_artists=[leg], bbox_inches='tight')
    plt.close()

    #RNAfold MFE vs Rollout MFE
    fig, ax = plt.subplots()
    line = np.linspace(-0.7, -0.1, 20)
    for fname in fnames:
        ax.scatter(all_data[RNAFOLD_MFE][families == fname] / all_data[SEQ_LEN][families == fname], all_data[ROLLOUT_MFE][families == fname] / all_data[SEQ_LEN][families == fname], label=fname, alpha=0.4)
    ax.scatter(all_data[RNAFOLD_MFE][best_rollout] / all_data[SEQ_LEN][best_rollout], all_data[ROLLOUT_MFE][best_rollout] / all_data[SEQ_LEN][best_rollout], c='red', s=8, label="Best ExpertRNA")
    ax.scatter(all_data[RNAFOLD_MFE][best_RNAfold] / all_data[SEQ_LEN][best_RNAfold], all_data[ROLLOUT_MFE][best_RNAfold] / all_data[SEQ_LEN][best_RNAfold], c='cyan', s=4, label="Best RNAfold")
    ax.plot(line, line, c='k', linewidth=0.75)
    ax.set_xlabel("RNAfold MFE (normalized)")
    ax.set_ylabel("ExpertRNA MFE (normalized)")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axis('scaled')
    leg = ax.legend(fontsize=14, bbox_to_anchor=(1.04, 1))
    fig.tight_layout()
    plt.savefig("RNAfold_roll_MFE_comparison_{}.pdf".format(out_number), bbox_extra_artists=[leg], bbox_inches='tight')
    plt.close()

        #how accurate is foldability??
    #plt.figure()
    #for fname in fnames:
    #    model = LinearRegression()
    #    model.fit(all_data[ROLLOUT_FOLDABILITY][families == fname].reshape(-1, 1), all_data[DATA_ROLL][families == fname].reshape(-1, 1))
    #    mi = min(all_data[ROLLOUT_FOLDABILITY])
    #    ma = max(all_data[ROLLOUT_FOLDABILITY])
    #    modelX =  np.linspace(mi, ma, 100)
    #    modelY = model.predict(modelX[:, np.newaxis])
    #    plt.scatter(all_data[ROLLOUT_FOLDABILITY][families == fname], all_data[DATA_ROLL][families == fname], label=fname+' {:.2f}'.format(model.score(all_data[ROLLOUT_FOLDABILITY][families == fname].reshape(-1, 1), all_data[DATA_ROLL][families == fname].reshape(-1, 1))))
    #    plt.plot(modelX, modelY)
    #plt.xlabel("ENTRNA foldability")
    #plt.ylabel("MCC")
    #plt.legend(fontsize=14)
    #plt.tight_layout()
    #plt.savefig("foldability_{}.pdf".format(out_number))
    #plt.close()

    #plt.figure(figsize=(len(all_data[DATA_FOLD])/18, 20))
    #ma = np.max(np.array([np.max(all_data[RNAFOLD_FOLDABILITY]), np.max(all_data[ROLLOUT_FOLDABILITY])]))
    #mi = np.min(np.array([np.min(all_data[RNAFOLD_FOLDABILITY]), np.min(all_data[ROLLOUT_FOLDABILITY])]))
    ##cmap = mpl.cm.viridis
    ##norm = mpl.colors.Normalize(vmin=mi, vmax=ma)
    #a = plt.scatter(range(len(all_data[DATA_FOLD])), all_data[DATA_FOLD], c=all_data[RNAFOLD_FOLDABILITY], cmap='viridis', vmin=mi, vmax=ma, marker='o', label="RNAfold")
    #b = plt.scatter(range(len(all_data[DATA_FOLD])), all_data[DATA_ROLL], c=all_data[ROLLOUT_FOLDABILITY], cmap='viridis', vmin=mi, vmax=ma, marker='v', label="Rollout")
    #plt.legend(handles=[a, b], fontsize=14)
    #plt.colorbar(label='Free energy')
    #plt.xlabel("Structure")
    #plt.ylabel("MCC")
    #plt.xticks(range(len(all_data[DATA_FOLD])), [n.replace('_', r'\_') for n in names], rotation='vertical', fontsize=6)
    #plt.xlim((-0.5, len(all_data[DATA_FOLD])+0.5))
    #plt.ylim((-1.1, 1.1))
    #plt.tight_layout() 
    #plt.savefig("all_data_{}.pdf".format(out_number))
    #plt.close()

if __name__ == "__main__":
    print('branch,avg_fold,avg_mcc,median_mcc,avg_improvement,stdev_improvement,pvaule')
    for i in range(5):
        main(i)
    #main(0)
