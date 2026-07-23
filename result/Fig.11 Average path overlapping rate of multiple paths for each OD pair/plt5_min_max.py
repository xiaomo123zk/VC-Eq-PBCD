import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import matplotlib.font_manager as fm

# 设置字体为新罗马字体
font_path = 'C:/Windows/Fonts/times.ttf'  # 替换为新罗马字体的路径
font_prop = fm.FontProperties(fname=font_path)

# 读取CSV文件，假设文件名为'data.csv'
data = pd.read_csv('Phi_ODin_new_PathOverlap4.csv', header=None)

# 绘制直方图
plt.hist(data.iloc[:, 0], bins=50, log=True) #,density=True,
plt.xlabel('OD pair paths overlap')
plt.ylabel('Cumulative OD pairs') 
plt.yscale('log')
# 格式化横轴的单位为科学计数
# plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
plt.tight_layout()
plt.savefig('Phi_OD_pathoverlap_Frequency2.png', dpi=300)  # 设置dpi为300，保存为高质量图片
plt.show()
