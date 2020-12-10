import pandas as pd
import matplotlib.pyplot as plt

def chi_plot(df, chi_name):
    cols = df.columns #['ILE182','MET170','PHE144','PHE105']
    N = len(cols)
    fig, axes = plt.subplots(nrows=5, ncols=1) 
    for i, col in enumerate(cols):
        print(i)
        df[col].plot(ax=axes[i], figsize=(20,40), fontsize=20)
        axes[i].set_ylim(-200,200)
        axes[i].set_title(label=f'{col}-{chi_name}', fontsize=20)
        fig.tight_layout()    

        if (i+1)%5 == 0: 
            plt.savefig(f'{i}.png')
            fig, axes = plt.subplots(nrows=5, ncols=1) #initialisation

def main():
    dfX1 = pd.read_csv('X1.csv', index_col = 0)
    dfX2 = pd.read_csv('X2.csv', index_col = 0)

    #cols = ['PHE12', 'PHE105', 'PHE123', 'TRP143','PHE144','PHE146','TYR173']
    chi_plot(dfX1, 'X1')
    #chi_plot(dfX2, 'X2')

    #df.plot.kde(figsize=(20,10), linewidth=5)
    #plt.legend(loc=2, prop={'size': 30})
    #plt.xlim(-180,180)
    plt.show()
if __name__ == "__main__":
    main()




