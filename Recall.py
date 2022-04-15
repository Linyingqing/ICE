# encoding=utf-8
import pandas as pd
data = pd.DataFrame(pd.read_csv('Prostate_Total_recall.csv'))
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams["font.family"] = 'Times New Roman'  #默认字体类型
mpl.rcParams["axes.linewidth"] = 1   #轴线边框粗细（默认的太粗了）
mpl.rcParams["font.size"] = 22

x = data['TOP']
ICE = data['ICE']
RLAG = data['RLAG']
EntroRank = data['EntroRank']
Dawnrank = data['Dawnrank']
Muf_max = data['Muf_max']
Muf_sum = data['Muf_sum']
Diffusion = data['Diffusion']
fig, ax = plt.subplots(1, 1,figsize=(8,7))

plt.ylim(0, 0.144)  # 限定横轴的范围 pro
#plt.ylim(0, 0.088)  # 限定横轴的范围bre
#plt.ylim(0, 0.090)  # 限定横轴的范围 lung
plt.xlim(0, 205)  # 限定纵轴的范围

plt.plot(x, Muf_sum, marker=',', linestyle='-',label='Muff_Sum',color='Purple')
plt.plot(x, Muf_max, marker=',', linestyle='-',label='Muff_Max',color='black')
plt.plot(x, Diffusion, marker=',', linestyle='-',label='Diffusion',color='Cyan')
plt.plot(x, Dawnrank, marker=',', linestyle='-',label='Dawnrank',color = 'forestgreen')
plt.plot(x, EntroRank, marker=',', linestyle='-',label='EntroRank',color = 'Orange')
plt.plot(x, RLAG, marker=',', linestyle='-',label='RLAG',color = 'blue')
plt.plot(x, ICE, marker=',', linestyle='-',label='ICE',color='r')
plt.legend(loc="lower right",fontsize=15)  # 让图例生效

plt.margins(0)
plt.subplots_adjust(bottom=0.15)
plt.xlabel('Top N Genes',fontsize=24)
plt.ylabel("Recall",fontsize=24)
plt.title("Prostate",fontsize=32,fontweight='bold',pad=20) #标题

axins = ax.inset_axes((0.12, 0.66, 0.4, 0.3))
axins.plot(x, Muf_sum, marker=',', linestyle='-',label='Muff_Sum',color='Purple')
axins.plot(x, Muf_max, marker=',', linestyle='-',label='Muff_Max',color='black')
axins.plot(x, Diffusion, marker=',', linestyle='-',label='Diffusion',color='Cyan')
axins.plot(x, Dawnrank, marker=',', linestyle='-',label='Dawnrank',color = 'forestgreen')
axins.plot(x, EntroRank, marker=',', linestyle='-',label='Muff_max',color='Orange')
axins.plot(x, RLAG, marker=',', linestyle='-',label='RLAG',color = 'blue')
axins.plot(x, ICE, marker=',', linestyle='-',label='ICE',color='r')

# 调整子坐标系的显示范围
axins.set_xlim(0, 31)

axins.set_ylim(0, 0.031)#pro
#axins.set_ylim(0, 0.0188)#lung
#axins.set_ylim(0, 0.0170)#bre
plt.tight_layout()
plt.savefig('figure/pro_rec.jpeg', dpi=900)
plt.savefig('figure/pro_rec.eps', dpi=900)
plt.savefig('figure/pro_rec.pdf', dpi=900)


plt.show()
