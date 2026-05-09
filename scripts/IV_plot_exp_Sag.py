import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from pathlib import Path

def plot_anode_iv_from_csv(
    csv_file: str,
    output_dir: str = "figure",
    basename: str = "anode_IV_Sag",
    dpi: int = 300,
    square_size: float = 3.3,    # 共有コードの雰囲気に合わせ小さめ正方形
    formats=("png", "pdf", "eps", "jpg"),
):
    """
    CSV仕様：
      1列目: 実験 Anode Current
      2列目: 実験 Anode Voltage
      3,4列目: 計算結果1 (Anode Current, Anode Voltage)
      5,6列目: 計算結果2 (Anode Current, Anode Voltage)
      ...
    ・1行目にヘッダあり
    ・凡例は偶数列（2,4,6,...列目）のヘッダ名を使用
    """

    # ===== スタイル設定（共有コードに準拠） =====
    # mathtext を Arial、軸ラベルは太字
    mpl.rcParams['mathtext.fontset'] = 'custom'
    mpl.rcParams['mathtext.rm'] = 'Arial'
    mpl.rcParams['mathtext.it'] = 'Arial:italic'
    mpl.rcParams['mathtext.bf'] = 'Arial:bold'
    mpl.rcParams['axes.labelweight'] = 'bold'

    # ===== データ読込 =====
    df = pd.read_csv(csv_file)
    ncols = df.shape[1]
    if ncols < 2:
        raise ValueError("CSVは少なくとも2列（Current, Voltage）が必要です。")

    # 出力先
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ===== 図作成（単一・正方形）=====
    fig, ax = plt.subplots(figsize=(square_size, square_size))

    # 軸ラベルと範囲（共有コードのサイズ感：fontsize=10, labelpad=2）
    ax.set_xlabel("Anode Current, A", fontsize=10, labelpad=2)
    ax.set_ylabel("Anode Voltage, V", fontsize=10, labelpad=2)
    ax.set_xlim(0.0, 0.4)   # ご指定の更新
    ax.set_ylim(20, 55)

    # 目盛：上下左右すべてに主・補助目盛、内向き。ラベルサイズは共有コード準拠で 8pt
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', which='major', labelsize=8, direction='in', top=True, right=True)
    ax.tick_params(axis='both', which='minor', direction='in', top=True, right=True)

    # 正方形アスペクト維持
    if hasattr(ax, "set_box_aspect"):
        ax.set_box_aspect(1)
    else:
        ax.set_aspect(1 / ax.get_data_ratio(), adjustable="box")

    # ===== 実験データ（1,2列目）：黒実線・マーカーなし。凡例ラベルは2列目ヘッダ =====
    x_exp = df.iloc[:, 0]
    y_exp = df.iloc[:, 1]
    label_exp = str(df.columns[1])  # 偶数列ヘッダ
    ax.plot(
        x_exp, y_exp,
        color="black", linestyle="-", linewidth=0.8, marker=None,
        label=label_exp
    )

    # ===== 計算結果（偶数列 4列目以降）：色付きマーカーのみ（線なし）=====
    # 先頭3本は red, blue, green、以降は C0, C1, ... を順に使用
    comp_color_seq = ["red", "C1", "blue"] + [f"C{i}" for i in range(10)]
    marker_list = ["o", "P", "s", "D", "v", "P", "X", "*", "<", ">", "h", "H",
                   "1", "2", "3", "4", "8", "p", "d"]

    comp_idx = 0
    for even_col_idx in range(3, ncols, 2):  # 0-based: 3,5,7,... (1-based: 4,6,8,...)
        xcol_idx = even_col_idx - 1
        if xcol_idx >= ncols:
            break

        x = df.iloc[:, xcol_idx]
        y = df.iloc[:, even_col_idx]
        label = str(df.columns[even_col_idx])

        color = comp_color_seq[comp_idx % len(comp_color_seq)]
        marker = marker_list[comp_idx % len(marker_list)]
        comp_idx += 1

        ax.plot(
            x, y,
            linestyle="None",                # 線なし
            marker=marker, markersize=4,     # 共有コードに合わせて小さめ
            markerfacecolor=color, markeredgecolor=color,
            label=label
        )

    # 凡例：共有コード準拠（fontsize=7, 枠なし）
    ax.legend(loc="lower right", frameon=True, fontsize=7)

    fig.tight_layout()

    # 保存
    for fmt in formats:
        fig.savefig(outdir / f"{basename}.{fmt}", dpi=dpi)

    plt.show()


if __name__ == "__main__":
    # 例: plot_anode_iv_from_csv("anode_data.csv")
    plot_anode_iv_from_csv("IV_plot_exp_Sag_rev.csv")
