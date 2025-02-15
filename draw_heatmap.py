import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.transforms import Bbox

def yelab_palette():
    # Define the color palette to match Gnuplot
    colors = [
        (1, 1, 1),      # White at value 0
        (0.33, 0, 0),   # Dark red at value 1
        (0, 1, 0)       # Green at value 6
    ]
    cmap1 = mcolors.LinearSegmentedColormap.from_list("custom_palette", colors[:2], N=256)
    cmap2 = mcolors.LinearSegmentedColormap.from_list("custom_palette", colors[1:], N=256)

    # """Combine two colormaps with a given split ratio.""" 1/6 customaized !
    colors1 = cmap1(np.linspace(0, 1, int(256 * 1/6)))  # Colors from the first colormap
    colors2 = cmap2(np.linspace(0, 1, 256 - len(colors1)))  # Colors from the second colormap
    combined_colors = np.vstack((colors1, colors2))  # Stack the two color arrays


    return mcolors.LinearSegmentedColormap.from_list("combined_cmap", combined_colors)

def extent2inch(extent, dpi, ax):
    """在plt.matplotlib.pyplot使用中，我如何得知ax.imshow（）中色块的以英尺为单位的bbox值"""
    left, right, top, bottom = extent
    x0_pixel, y0_pixel = ax.transData.transform((left, bottom))  # 左下角像素坐标
    x1_pixel, y1_pixel = ax.transData.transform((right, top))  # 右上角像素坐标

    # 计算像素宽度和高度
    width_pixel = x1_pixel - x0_pixel
    height_pixel = y1_pixel - y0_pixel

    # 将像素转换为英寸
    width_inch = width_pixel / dpi
    height_inch = height_pixel / dpi

    return width_inch*10 + 0.3, height_inch*10 + 0.2


def heatmap(file_path):
    # Load data from the file (replace 'ribo_mass.txt' with your file path)

    data = pd.read_csv(file_path, sep='\t')
    name= file_path[:-4]
    data['NAF'] = data['NAF'].astype(float)
    data['cheat'] = np.zeros((len(data['NAF']),1))
    # print(data['NAF'])


    # Generate the combined colormap
    combined_cmap = yelab_palette()

    # Normalize values to match the color bar range [0:6]
    norm = mcolors.Normalize(vmin=0, vmax=6)

    # Create the figure and axis
    scale=1
    fig, ax = plt.subplots(figsize=(4*scale, 8*scale))
    ax.axis('off')
    # Plot the heatmap
    heatdata = data[['NAF','cheat']]
    heatmap = ax.imshow(heatdata, cmap=combined_cmap, norm=norm, aspect=0.64, origin="upper")

    for i in range(heatdata.shape[0]):  # Iterate over rows
        for j in range(heatdata.shape[1]):  # Iterate over columns
            ax.text(0.64*scale, i*scale, f"{data['Protein'][i]}",  # Text to display
                    ha="left", va="center",  # Right alignment
                    color="black",
                    fontfamily='Arial',
                    fontsize=3)  # Adjust color for contrast

    # cluster name
    ax.text(0.5, heatdata.shape[0]+0.5, f"{name}",  # Text to display
            ha="center", va="center",  # Right alignment
            color="black",
            fontfamily='Arial',
            fontsize=3)  # Adjust color for contrast

    # Add a colorbar
    # cbar = plt.colorbar(heatmap, ax=ax, orientation="vertical")
    # cbar.set_label("Intensity")

    # Set the axis ranges (xrange and yrange)
    ax.set_xlim(-0.5, 35*scale)  # xrange [-0.5:35]
    ax.set_ylim(-0.5, 120*scale)  # yrange [-0.5:120]

    # Set axis labels
    # ax.set_xlabel("X-axis (Columns)")
    # ax.set_ylabel("Y-axis (Rows)")

    # Set the title
    # ax.set_title(name)

    # output
    extent = heatmap.get_extent()
    x_width, y_width = extent2inch(extent, 800, ax)
    # print(x_width, y_width)

    # Save the plot to a png file & remove most blank region
    plt.tight_layout()
    bbox = Bbox([[0.1,0.1],[0.1 + x_width,0.1 + y_width]]) # (left, bottom, right, top)
    plt.savefig(f"{name}.png",dpi=800, bbox_inches=bbox)
    # plt.savefig(f"test.png", dpi=800)
    plt.close()

    print(f"Heatmap saved as '{name}.png'")
    # print(f"Heatmap saved as 'test.png'")


if __name__ == '__main__':
    heatmap('assembly_factors.txt')

    # 获取当前目录下所有文件名
    all_files = os.listdir()

    # 筛选出所有 .txt 文件
    # txt_files = [f for f in all_files if f.endswith('.txt')]
    # print(txt_files)
    #
    # for f in txt_files:
    #     try:
    #         heatmap(f)
    #     except:
    #         print(f"Error occurs in file: {f}...")

