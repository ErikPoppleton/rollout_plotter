import numpy as np
from re import finditer
from os import listdir
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.linear_model import LinearRegression

def calc_MCC(prediction, actual, name, alg):
    TP = 0
    TN = 0
    FP = 0
    FN = 0
    with open("{}_{}_colormap.col".format(name, alg), 'w+') as f:

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

def make_fasta_db(name, alg, seq, struct):
    with open("{}_{}.fasta".format(name, alg), 'w+') as f:
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
    names = []
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
                    actual = l[2]
                    RNAfold = l[3]
                    rollout = l[4]
                    family = name.split("_")[0]
                    actual_energy.append(float(l[5]))
                    RNAfold_energy.append(float(l[6]))
                    rollout_energy.append(float(l[7]))
                    
                    make_fasta_db(name, "actual", seq, actual)
                    make_fasta_db(name, "RNAfold", seq, RNAfold)
                    make_fasta_db(name, "rollout", seq, rollout)
                    actual = parse_dot_bracket(actual)            
                    rollout = parse_dot_bracket(rollout)
                    RNAfold = parse_dot_bracket(RNAfold)
                    acc_fold = calc_MCC(RNAfold, actual, name, "RNAfold")
                    acc_roll = calc_MCC(rollout, actual, name, "rollout")
                    data_fold.append(acc_fold)
                    data_roll.append(acc_roll)
                    data_diff.append(acc_roll - acc_fold)
                    names.append(name)
                    families.append(family)
                    #print(name, acc_fold, acc_roll, acc_roll-acc_fold)

    all_data = np.array([data_fold, data_roll, data_diff, actual_energy, RNAfold_energy, rollout_energy], dtype=float)
    DATA_FOLD = 0
    DATA_ROLL = 1
    DATA_DIFF = 2
    ACTUAL_ENERGY = 3    
    RNAFOLD_ENERGY = 4
    ROLLOUT_ENERGY = 5
    sort = all_data[DATA_DIFF].argsort()
    all_data = all_data[:,sort]
    names = np.array(names)[sort]
    families = np.array(families)[sort]


    print("{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}".format(out_number, np.mean(all_data[ROLLOUT_ENERGY]), np.mean(all_data[DATA_ROLL]), np.median(all_data[DATA_ROLL]), np.mean(all_data[DATA_DIFF])))

    with open("scores_{}.dat".format(out_number), 'w') as f:
        for n, d in zip(names, all_data.T):
            f.write("{} {} {} {}\n".format(n, d[DATA_FOLD], d[DATA_ROLL], d[DATA_DIFF]))


    mpl.rcParams.update({'font.size': 22})
    mpl.rc('text', usetex=True)

    #how accurate is foldability??
    plt.clf()
    plt.figure()
    fnames = sorted(list(set(families)))
    for fname in fnames:
        model = LinearRegression()
        model.fit(all_data[ROLLOUT_ENERGY][families == fname].reshape(-1, 1), all_data[DATA_ROLL][families == fname].reshape(-1, 1))
        mi = min(all_data[ROLLOUT_ENERGY])
        ma = max(all_data[ROLLOUT_ENERGY])
        modelX =  np.linspace(mi, ma, 100)
        modelY = model.predict(modelX[:, np.newaxis])
        plt.scatter(all_data[ROLLOUT_ENERGY][families == fname], all_data[DATA_ROLL][families == fname], label=fname+' {:.2f}'.format(model.score(all_data[ROLLOUT_ENERGY][families == fname].reshape(-1, 1), all_data[DATA_ROLL][families == fname].reshape(-1, 1))))
        plt.plot(modelX, modelY)
    plt.xlabel("ENTRNA foldability")
    plt.ylabel("MCC")
    plt.legend()
    plt.tight_layout()
    plt.savefig("foldability_{}.pdf".format(out_number))

    #comparison scatterplot
    plt.clf()
    plt.figure()
    for fname in fnames:
        plt.scatter(all_data[DATA_FOLD][families == fname], all_data[DATA_ROLL][families == fname], label=fname)
        plt.xlabel("RNAfold MCC")
        plt.ylabel("Rollout MCC")
        plt.legend()
        plt.tight_layout()
        plt.savefig("scatter_comparison_{}.pdf".format(out_number))

    #histograms!
    plt.clf()
    plt.figure()
    #overlay
    plt.hist(data_fold, bins=50, range=(-1, 1), alpha=0.4, label="RNAfold")
    plt.hist(data_roll, bins=50, range=(-1, 1), alpha=0.4, label="Rollout")
    plt.xlabel("Matthews Correlation Coefficient")
    plt.ylabel("Number")
    plt.legend()
    plt.tight_layout()
    plt.savefig("comparison_{}.pdf".format(out_number))

    #difference
    plt.clf()
    plt.hist(data_diff, bins=20)
    plt.xlabel("MCC difference (Rollout-RNAfold)")
    plt.ylabel("Number")
    plt.savefig("diff_{}.pdf".format(out_number))



    plt.clf()
    plt.figure(figsize=(200, 20))
    ma = np.max(np.array([np.max(all_data[RNAFOLD_ENERGY]), np.max(all_data[ROLLOUT_ENERGY])]))
    mi = np.min(np.array([np.min(all_data[RNAFOLD_ENERGY]), np.min(all_data[ROLLOUT_ENERGY])]))
    #cmap = mpl.cm.viridis
    #norm = mpl.colors.Normalize(vmin=mi, vmax=ma)
    a = plt.scatter(range(len(all_data[DATA_FOLD])), all_data[DATA_FOLD], c=all_data[RNAFOLD_ENERGY], cmap='viridis', vmin=mi, vmax=ma, marker='o', label="RNAfold")
    b = plt.scatter(range(len(all_data[DATA_FOLD])), all_data[DATA_ROLL], c=all_data[ROLLOUT_ENERGY], cmap='viridis', vmin=mi, vmax=ma, marker='v', label="Rollout")
    plt.legend(handles=[a, b])
    plt.colorbar(label='Free energy')
    plt.xlabel("Structure")
    plt.ylabel("MCC")
    plt.xticks(range(len(all_data[DATA_FOLD])), names, rotation='vertical', fontsize=6)
    plt.xlim((-0.5, len(all_data[DATA_FOLD])+0.5))
    plt.ylim((-1.1, 1.1))
    plt.tight_layout() 
    plt.tight_layout()   
    plt.savefig("all_data_{}.pdf".format(out_number))

if __name__ == "__main__":
    print('\tavg_fold\tavg_mcc\tmedian_mcc\tavg_improvement')
    #for i in range(5):
    #    main(i)
    main(1)
