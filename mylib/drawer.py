import matplotlib.pyplot as plt
from . import type

def d(data: type.OutputData, path):
    x = [p.x + p.ax for p in data.points]
    y = [p.y + p.ay for p in data.points]
    fig, ax = plt.subplots(figsize=(8, 10), tight_layout = True)
    ax.set_aspect('equal')
    
    scatter = ax.scatter(
    x,  # X坐标
    y,  # Y坐标
    c='blue',           # 点颜色
    s=40,               # 点大小
    alpha=0.7,          # 透明度
    edgecolors='black', # 边框颜色
    marker='o'          # 点形状
    )
    ax.plot(
        x,  # X坐标
        y,  # Y坐标
        color='blue',      # 线颜色
        linewidth=1.5,     # 线宽
        linestyle='-',     # 线型
    )
    x_min = -5#min(x) - 0.2
    y_min = -0.3#min(y) - 0.2
    x_max = 10.3#max(x) + 0.2
    y_max = 10#max(y) + 0.2
    ax.set_xlim(x_min,x_max)
    ax.set_ylim(y_min,y_max)
    # plt.show()
    fig.savefig(path)
    plt.close()