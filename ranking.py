import numpy as np
import matplotlib.pyplot as plt

def get_rankings_for(start, stop):
    s = slice(start, stop)
    names = []

    with open('scores_1.dat', 'r') as f:
        data = f.readlines()
        ranks = np.zeros((len(range(*s.indices(len(data)))),4))
        
        for i, line in enumerate(data[s]):
            line = line.split()
            names.append(line[0])
            if start >= 0:
                ranks[i][0] = start + i
            else:
                ranks[i][0] = len(data) + start + i

    for i in [2, 3, 4]:
        with open('scores_{}.dat'.format(i), 'r') as f:
            data = f.readlines()
            for j, line in enumerate(data):
                line = line.split()
                for k, name in enumerate(names):
                    if line[0] == name:
                        ranks[k][i-1] = j

    return(names, ranks)

if __name__ == "__main__":
    bad_names, bad_ranks = get_rankings_for(0, 10)
    good_names, good_ranks = get_rankings_for(-10, None)

    
    for n, r in zip(bad_names, bad_ranks):
        plt.plot([1, 2, 3, 4], r, label=n, marker='v')
    for n, r in zip(good_names, good_ranks):
        plt.plot([1, 2, 3, 4], r, label=n, marker='o')   
    plt.xlabel('Branch')
    plt.ylabel('MCC Rank')
    plt.xticks([1, 2, 3, 4])
    plt.legend(bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    plt.savefig('rank.png')

    all_names, all_ranks = get_rankings_for(0, None)
    plt.clf()
    plt.figure()
    plt.plot(np.arange(len(all_ranks)), np.std(all_ranks, axis=1))
    plt.xlabel("Initial Rank")
    plt.ylabel("Rank StDev")
    plt.tight_layout()
    plt.savefig('rank_dev.png')
    