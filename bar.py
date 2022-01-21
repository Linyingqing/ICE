
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams["font.family"] = 'Times New Roman'  
fig = plt.figure()
patterns = ["/ /"  , "\ \ " , "XXX","...","**" ]
ax = fig.add_subplot(111)
plt.ylim(ymax=150, ymin=0)

labels = ['Breast', 'Prostate', 'Lung']
TOP50_CERW = [41 , 42, 49]
TOP50_EntroRank = [37 , 37, 39]
TOP50_Dawnrank = [41 , 37, 35]
TOP50_Muff_Max = [28 , 23, 21]
TOP50_Muff_Sum = [29 , 34, 41]

TOP100_CERW = [77 , 72, 83]
TOP100_EntroRank = [68 , 62, 73]
TOP100_Dawnrank = [71 , 62, 67]
TOP100_Muff_Max = [57 , 42, 35]
TOP100_Muff_Sum = [62 , 60, 70]

TOP200_CERW = [140 , 117, 143]
TOP200_EntroRank = [125 , 110, 130]
TOP200_Dawnrank = [122 , 111, 121]
TOP200_Muff_Max = [90 , 89, 91]
TOP200_Muff_Sum = [118 , 117, 127]

index = np.arange(len(labels))
bar_width = 0.16

ax.bar(index,TOP200_CERW, bar_width,color='white', edgecolor='m', hatch=patterns[0], label='Top200_ICE')
ax.bar(index+bar_width,TOP200_EntroRank, bar_width,color='white', edgecolor='m', hatch=patterns[1], label='Top200_EntroRank')
ax.bar(index+2*bar_width,TOP200_Dawnrank, bar_width,color='white', edgecolor='m', hatch=patterns[2], label='Top200_Dawnrank')
ax.bar(index+3*bar_width,TOP200_Muff_Max, bar_width,color='white', edgecolor='m', hatch=patterns[3], label='Top200_Muff_Max')
ax.bar(index+4*bar_width,TOP200_Muff_Sum, bar_width,color='white', edgecolor='m', hatch=patterns[4], label='Top200_Muff_Sum')

ax.bar(index,TOP100_CERW, bar_width,color='white', edgecolor='b', hatch=patterns[0], label='TOP100_ICE')
ax.bar(index+bar_width,TOP100_EntroRank, bar_width,color='white', edgecolor='b', hatch=patterns[1], label='Top100_EntroRank')
ax.bar(index+2*bar_width,TOP100_Dawnrank, bar_width,color='white', edgecolor='b', hatch=patterns[2], label='Top100_Dawnrank')
ax.bar(index+3*bar_width,TOP100_Muff_Max, bar_width,color='white', edgecolor='b', hatch=patterns[3], label='Top100_Muff_Max')
ax.bar(index+4*bar_width,TOP100_Muff_Sum, bar_width,color='white', edgecolor='b', hatch=patterns[4], label='Top100_Muff_Sum')

ax.bar(index,TOP50_CERW, bar_width,color='white', edgecolor='r', hatch=patterns[0], label='TOP50_ICE')
ax.bar(index+bar_width,TOP50_EntroRank, bar_width,color='white', edgecolor='r', hatch=patterns[1], label='Top50_EntroRank')
ax.bar(index+2*bar_width,TOP50_Dawnrank, bar_width,color='white', edgecolor='r', hatch=patterns[2], label='Top50_Dawnrank')
ax.bar(index+3*bar_width,TOP50_Muff_Max, bar_width,color='white', edgecolor='r', hatch=patterns[3], label='Top50_Muff_Max')
ax.bar(index+4*bar_width,TOP50_Muff_Sum, bar_width,color='white', edgecolor='r', hatch=patterns[4], label='Top50_Muff_Sum')

box = ax.get_position()

ax.set_position([box.x0, box.y0, box.width , box.height* 1])
ax.legend(loc='center', bbox_to_anchor=(0.498, 1.2),ncol=3,fontsize=10.88)

plt.xticks(index+2*bar_width, labels,fontsize=16)
plt.yticks(np.arange(0, 150, 40),fontsize=12)

plt.tight_layout()

plt.savefig('pdf/bar.jpeg', dpi=900)
plt.savefig('pdf/bar.eps', dpi=900)
plt.savefig('pdf/bar.pdf', dpi=900)
plt.show()
