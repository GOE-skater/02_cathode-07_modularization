import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.colors import ListedColormap, LogNorm
import os
import sys

def find_column_index_by_header(filename: str, header_name: str) -> int:
    """CSVの1行目ヘッダーから header_name と完全一致（前後空白は無視）する列のindexを返す"""
    df0 = pd.read_csv(filename, nrows=0)  # ヘッダーだけ読む
    cols = [str(c).strip() for c in df0.columns.tolist()]
    target = header_name.strip()
    try:
        idx = cols.index(target)
        print(f"[info] header match: '{target}' -> column index {idx}")
        return idx
    except ValueError:
        raise ValueError(
            f"指定ヘッダー '{header_name}' が見つかりません。利用可能なヘッダー: {cols}"
        )

def format_tick_labels(vals, icon_exp: int, ndigits: int):
    """カラーバー用の数値ラベルを icon_exp/ndigits 設定でフォーマット"""
    if icon_exp == 1:
        fmt = f"{{:.{ndigits}e}}"   # 指数表記
    else:
        fmt = f"{{:.{ndigits}f}}"   # 少数表記
    return [fmt.format(float(v)) for v in vals]

def plot_combined_streamlines_single_no_streamline(file_path):
    # Parameters
    filename = file_path
    xCol = 2              # x の列番号（0始まり）
    yCol = 3              # y の列番号（0始まり）

    # ▼これは上書きされる（output_name と一致する列を自動検出して dataCol に代入）
    dataCol = None

    # ▼デフォルトの手動レンジ（icon_vauto=0のときに使用）
    vmin = 0.0
    vmax = 45

    icon_vauto = 1        # 1: 値の範囲を自動設定
    icon_log = 0          # 1: カラーレンジを対数スケールに切替

    # ▼新規追加（カラーバーの数字表記）
    icon_exp = 0          # 1: 指数表記, 0: 少数表記
    ndigits = 1           # 少数点以下の桁数（デフォルト1）

    num_colors = 10       # 離散的な色数
    num_cbar_label = 5    # カラーバーのラベルの数
    cbar_label = r'$Plasma\ Potential\ (V)$'
    output_name = 'phi'  # 出力ファイル名 兼 「抽出したいヘッダー名」

    # === ここで output_name と同じヘッダーを持つ列番号を自動取得 ===
    dataCol = find_column_index_by_header(filename, output_name)

    xRange = [0, 0.0388]
    yRange = [0.00, 0.014]
    offset = 7.8e-3
    MaskRange = [
        [0.00, 10e-3 + offset, 0.000, 0.8e-3],
        [0.00, 5.5e-3 + offset, 0.8e-3, 4e-3],
        [0.00, offset, 4e-3, 9e-3],
        [0.00, 15e-3 + offset, 9e-3, 14e-3],
        [15e-3 + offset, 20e-3 + offset, 2e-3, 14e-3]
    ]

    gridDensity = 1000
    majorTickInterval = 0.004

    plt.rcParams["font.size"] = 8
    fig, ax = plt.subplots(figsize=(6, 3))

    im, vmin_used, vmax_used, is_log = plot_data(
        ax, filename, xCol, yCol, dataCol,
        xRange, yRange, MaskRange,
        vmin, vmax, icon_vauto, icon_log,
        gridDensity, majorTickInterval,
        num_colors, output_name
    )

    # --- カラーバー ---
    cbar = plt.colorbar(im, ax=ax, orientation='vertical', aspect=6, shrink=0.6)
    cbar.set_label(cbar_label, fontsize=10)

    # 線形/対数で等間隔tickを作成し、表記を icon_exp / ndigits で切替
    if is_log:
        ticks = np.logspace(np.log10(vmin_used), np.log10(vmax_used), num_cbar_label)
    else:
        ticks = np.linspace(vmin_used, vmax_used, num_cbar_label)

    cbar.set_ticks(ticks)
    cbar.set_ticklabels(format_tick_labels(ticks, icon_exp, ndigits))

    print("Colorbar Range (used) = ", vmin_used, vmax_used)

    plot_dir = 'figure'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    output_file = os.path.join(plot_dir, output_name)
    plt.savefig(f'{output_file}.png', dpi=300)
    plt.savefig(f'{output_file}.jpg', dpi=300)
    plt.savefig(f'{output_file}.pdf')
    plt.savefig(f'{output_file}.eps')

def plot_data(
    ax, filename, xCol, yCol, dataCol,
    xRange, yRange, MaskRange,
    vmin, vmax, icon_vauto, icon_log,
    gridDensity, majorTickInterval,
    num_colors, output_name
):
    # ▼ヘッダーありで全体を読む（1行目が列名）
    data = pd.read_csv(filename, header=0)

    x = data.iloc[:, xCol].to_numpy()
    y = data.iloc[:, yCol].to_numpy()
    values = data.iloc[:, dataCol].to_numpy()

    # グリッド化
    X, Y = np.meshgrid(
        np.linspace(xRange[0], xRange[1], gridDensity),
        np.linspace(yRange[0], yRange[1], gridDensity)
    )
    Z = griddata((x, y), values, (X, Y), method='linear')

    # マスク
    mask = np.zeros_like(Z, dtype=bool)
    for x_min, x_max, y_min, y_max in MaskRange:
        mask |= (X >= x_min) & (X <= x_max) & (Y >= y_min) & (Y <= y_max)
    Z[mask] = np.nan

    # 自動レンジ
    if icon_vauto == 1:
        if icon_log == 1:
            pos = Z[(Z > 0) & np.isfinite(Z)]
            if pos.size == 0:
                raise ValueError("ログスケールにできる正の値が見つかりません。icon_log=0 で実行してください。")
            #vmin = np.nanmin(pos)
            vmax = np.nanmax(pos)
        else:
            #vmin = np.nanmin(Z)
            vmax = np.nanmax(Z)
        print("Value Range (auto) = ", vmin, vmax)
    else:
        if icon_log == 1:
            if not (vmin is not None and vmin > 0):
                pos_min = np.nanmin(Z[(Z > 0) & np.isfinite(Z)])
                vmin = pos_min
            if vmax is None or vmax <= vmin:
                vmax = np.nanmax(Z[(Z > 0) & np.isfinite(Z)])

    # カラーマップ（NaNは白）
    cmap = plt.cm.get_cmap('jet', num_colors)
    cmap_with_white = cmap(np.linspace(0, 1, cmap.N))
    cmap_with_white = ListedColormap(cmap_with_white)
    cmap_with_white.set_bad(color='white')

    if icon_log == 1:
        vmin_safe = max(np.nextafter(0, 1), float(vmin))
        norm = LogNorm(vmin=vmin_safe, vmax=float(vmax))
        im = ax.imshow(
            Z, extent=xRange + yRange, origin='lower',
            alpha=1.0, cmap=cmap_with_white, norm=norm
        )
    else:
        im = ax.imshow(
            Z, extent=xRange + yRange, origin='lower',
            alpha=1.0, cmap=cmap_with_white,
            vmin=float(vmin), vmax=float(vmax)
        )

    ax.set_xlim(xRange)
    ax.set_ylim(yRange)
    ax.set_xlabel('axial position (m)', fontsize=10, labelpad=5)
    ax.set_ylabel('radial position (m)', fontsize=10, labelpad=5)

    ax.set_xticks(np.arange(xRange[0], xRange[1] + 1e-3, majorTickInterval))
    ax.set_yticks(np.arange(yRange[0], yRange[1] + 1e-3, majorTickInterval))
    ax.tick_params(axis='both', which='both', direction='out')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')

    ax.axis('tight')
    ax.set_aspect('equal', adjustable='box')

    return im, float(vmin), float(vmax), int(icon_log) == 1

if __name__ == '__main__':
    args = sys.argv
    file_path = args[1]
    plot_combined_streamlines_single_no_streamline(file_path)
