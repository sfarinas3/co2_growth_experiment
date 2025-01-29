import matplotlib.pyplot as plt
import seaborn as sn

def plot_boxplots(df, x, y, hue, title):
    ax = sn.boxplot(data = df, x=x, y=y, hue=hue)
    plt.title(title)
    plt.show()

def plot_multiple_boxplots(df, by, x, y, hue, subplots = {'nrows':1, 'ncols':2, 'figsize':(None,None)}):
    df = df.dropna(subset=[y])
    by_labels = df[by].unique()
    fig, axes = plt.subplots(**subplots)
    cnt = 0
    while cnt < len(by_labels):
        for irow in range(subplots['nrows']):
            for icol in range(subplots['ncols']):
                label = by_labels[cnt]
                dftmp = df[df[by] == label]
                axi = axes[irow,icol]
                sn.boxplot(data = dftmp, x=x, y=y, hue=hue, ax = axi)
                axi.set_title(label, weight = 'bold', fontsize = 13)
                axi.set_xlabel('')
                axi.set_ylabel('')
                cnt = cnt+1
    fig.text(0.5, 0.05, x, ha='center', fontsize = 12, weight = 'bold')
    fig.text(0.07, 0.5, y, va='center', rotation='vertical', fontsize = 12, weight = 'bold')
    plt.show()

def plot_sidebyside_boxplots(df, by, x, y1, y2, hue):
    df = df.dropna(subset=[y1])
    df = df.dropna(subset=[y2])
    by_labels = df[by].unique()
    for label in by_labels: 
        dftmp = df[df[by] == label]
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
        bp1 = sn.boxplot(data = dftmp, x=x, y=y1, hue=hue, ax = axes[0])
        bp2 = sn.boxplot(data = dftmp, x=x, y=y2, hue=hue, ax = axes[1])
        bp1.legend([],[], frameon=False)
        for ax in axes:
            ax.set_xlabel('')
        fig.text(0.5, 0.01, x, ha='center', fontsize = 12)
        fig.suptitle(label, fontsize = 13, weight = 'bold')
    plt.show()