import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.patches import Rectangle
import os
import sys

def find_column_index_by_header(filename: str, header_name: str) -> int:
    """CSVの1行目ヘッダーから header_name と完全一致（前後空白は無視）する列のindexを返す"""
    df0 = pd.read_csv(filename, nrows=0)  # ヘッダーだけ読む
    cols = [str(c).strip() for c in df0.columns.tolist()]
    target = str(header_name).strip()
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

def draw_mask_regions(ax, MaskRange, zorder=3):
    """マスク領域を白で塗りつぶす"""
    for x_min, x_max, y_min, y_max in MaskRange:
        ax.add_patch(
            Rectangle(
                (x_min, y_min),
                x_max - x_min,
                y_max - y_min,
                facecolor='white',
                edgecolor='none',
                zorder=zorder
            )
        )

def plot_combined_streamlines_single_no_streamline(file_path):
    # Parameters
    filename = file_path
    xCol = 2              # x の列番号（0始まり）
    yCol = 3              # y の列番号（0始まり）

    # ▼コンター値はCSV列ではなくベクトルのノルムから計算するため dataCol は使用しない
    dataCol = None

    # ▼デフォルトの手動レンジ（icon_vauto=0のときに使用）
    vmin = 1e1
    vmax = 1e4

    icon_vauto = 0        # 1: 値の範囲を自動設定
    icon_log = 1          # 1: カラーレンジを対数スケールに切替
    icon_contour = 1      # 1: コンター表示, 0: ベクトルのみ（カラーバー無し）

    # ▼カラーバー書式
    icon_exp = 1          # 1: 指数表記, 0: 少数表記
    ndigits = 1           # 少数点以下の桁数（デフォルト1）

    # ▼ベクトル表示
    icon_vector = 2       # 2: 流線図 1: ベクトルを重ねて描画, 0: 描画しない
    vector_size = 0.005   # 単位ベクトルに対して掛けるスケール
    # ランダム間引き設定
    vector_keep_ratio = 0.0015
    vector_max_points = 1500
    vector_random_seed = 42

    # ▼CSVヘッダー名（必要に応じて変更してください）
    vec_x_header = 'Uix'    # x方向成分
    vec_y_header = 'Uir'    # y方向成分
    vec_z_header = 'Uip'    # ★ z方向成分（任意の列名を指定可能）

    num_colors = 10
    num_cbar_label = 5
    cbar_label = r'$Ion\ Velocity\ (m/s)$'
    output_name = 'Ui_mag'  # 出力ファイル名（列は使わないが名前はそのまま）

    # === ベクトル列インデックス（★3成分すべて取得） ===
    vec_x_col = find_column_index_by_header(filename, vec_x_header)
    vec_y_col = find_column_index_by_header(filename, vec_y_header)
    vec_z_col = find_column_index_by_header(filename, vec_z_header)

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
    ax.set_facecolor('#f0f0f0')

    im, vmin_used, vmax_used, is_log = plot_data(
        ax, filename, xCol, yCol, dataCol,  # dataColは未使用
        xRange, yRange, MaskRange,
        vmin, vmax, icon_vauto, icon_log, icon_contour,
        gridDensity, majorTickInterval,
        num_colors, output_name,
        # ▼ベクトル描画・ノルム作成に使用
        icon_vector, vec_x_col, vec_y_col, vec_z_col, vector_size,
        vector_keep_ratio, vector_max_points, vector_random_seed
    )

    # --- カラーバー（コンター表示時のみ） ---
    if im is not None and icon_contour == 1:
        cbar = plt.colorbar(im, ax=ax, orientation='vertical', aspect=6, shrink=0.6)
        cbar.set_label(cbar_label, fontsize=10)

        if is_log:
            ticks = np.logspace(np.log10(vmin_used), np.log10(vmax_used), num_cbar_label)
        else:
            ticks = np.linspace(vmin_used, vmax_used, num_cbar_label)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(format_tick_labels(ticks, icon_exp, ndigits))
        print("Colorbar Range (used) = ", vmin_used, vmax_used)
    else:
        print("Contour disabled: colorbar is not shown.")

    plot_dir = 'figure'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    output_file = os.path.join(plot_dir, output_name)
    plt.savefig(f'{output_file}.png', dpi=300)
    plt.savefig(f'{output_file}.jpg', dpi=300)
    plt.savefig(f'{output_file}.pdf')
    plt.savefig(f'{output_file}.eps')

def plot_data(
    ax, filename, xCol, yCol, dataCol,  # dataColは未使用（互換のため引数は残す）
    xRange, yRange, MaskRange,
    vmin, vmax, icon_vauto, icon_log, icon_contour,
    gridDensity, majorTickInterval,
    num_colors, output_name,
    # ▼ベクトル描画・ノルム作成に使用
    icon_vector, vec_x_col, vec_y_col, vec_z_col, vector_size,
    vector_keep_ratio, vector_max_points, vector_random_seed
):
    # ▼ヘッダーありで全体を読む（1行目が列名）
    data = pd.read_csv(filename, header=0)

    # 座標
    x = data.iloc[:, xCol].to_numpy()
    y = data.iloc[:, yCol].to_numpy()

    # --- ベクトル成分を取得 ---
    u_vals = data.iloc[:, vec_x_col].to_numpy()  # Ux
    v_vals = data.iloc[:, vec_y_col].to_numpy()  # Uy
    w_vals = data.iloc[:, vec_z_col].to_numpy()  # Uz（★新規）

    # --- 背景用のグリッドを作成 ---
    X, Y = np.meshgrid(
        np.linspace(xRange[0], xRange[1], gridDensity),
        np.linspace(yRange[0], yRange[1], gridDensity)
    )

    # --- ベクトルをグリッドへ補間 ---
    U = griddata((x, y), u_vals, (X, Y), method='linear')
    V = griddata((x, y), v_vals, (X, Y), method='linear')
    W = griddata((x, y), w_vals, (X, Y), method='linear')

    # --- ノルムを3成分で計算（★変更点：sqrt(U^2 + V^2 + W^2)） ---
    Z = np.sqrt(U**2 + V**2 + W**2)

    # マスク適用（Z, U, V, Wすべて）
    mask = np.zeros_like(Z, dtype=bool)
    for x_min, x_max, y_min, y_max in MaskRange:
        mask |= (X >= x_min) & (X <= x_max) & (Y >= y_min) & (Y <= y_max)
    for arr in (Z, U, V, W):
        arr[mask] = np.nan

    # 自動レンジ（コンター表示時のみ必要だが、値は返しておく）
    if icon_vauto == 1:
        if icon_log == 1:
            pos = Z[(Z > 0) & np.isfinite(Z)]
            if pos.size == 0 and icon_contour == 1:
                raise ValueError("ログスケールにできる正の値が見つかりません。icon_log=0 で実行してください。")
            if pos.size > 0:
                vmin = float(np.nanmin(pos))
                vmax = float(np.nanmax(pos))
        else:
            vmin = float(np.nanmin(Z))
            vmax = float(np.nanmax(Z))
        print("Value Range (auto) = ", vmin, vmax)
    else:
        if icon_log == 1:
            if not (vmin is not None and vmin > 0):
                pos_min = np.nanmin(Z[(Z > 0) & np.isfinite(Z)])
                vmin = float(pos_min)
            if vmax is None or vmax <= vmin:
                vmax = float(np.nanmax(Z[(Z > 0) & np.isfinite(Z)]))

    # カラーマップ（NaNは白）
    cmap = plt.cm.get_cmap('jet', num_colors)
    cmap_with_white = cmap(np.linspace(0, 1, cmap.N))
    cmap_with_white = ListedColormap(cmap_with_white)
    cmap_with_white.set_bad(color='white')  # NaN（マスク）は白

    im = None
    # 背景コンター描画（必要なときだけ）
    if icon_contour == 1:
        if icon_log == 1:
            vmin_safe = max(np.nextafter(0, 1), float(vmin))
            norm = LogNorm(vmin=vmin_safe, vmax=float(vmax))
            im = ax.imshow(
                Z, extent=xRange + yRange, origin='lower',
                alpha=1.0, cmap=cmap_with_white, norm=norm, zorder=1
            )
        else:
            im = ax.imshow(
                Z, extent=xRange + yRange, origin='lower',
                alpha=1.0, cmap=cmap_with_white,
                vmin=float(vmin), vmax=float(vmax), zorder=1
            )

    # --- マスク領域を白で塗り重ね（コンター有無に関わらず） ---
    draw_mask_regions(ax, MaskRange, zorder=3)

    # --- ベクトル重ね描画（★x–y 成分のみを使用） ---
    if icon_vector == 1:
        # マスク・NaN除外
        good = np.isfinite(U) & np.isfinite(V)
        idx_all = np.flatnonzero(good)
        if idx_all.size > 0:
            keep_ratio = float(np.clip(vector_keep_ratio, 0.0, 1.0))
            if keep_ratio <= 0.0:
                n_keep = 0
            else:
                n_keep = int(np.ceil(idx_all.size * keep_ratio))
                if vector_max_points is not None:
                    n_keep = min(n_keep, int(vector_max_points))
                n_keep = max(n_keep, 1)

            if n_keep > idx_all.size:
                n_keep = idx_all.size

            rng = np.random.default_rng(vector_random_seed) if vector_random_seed is not None else np.random.default_rng()
            chosen = rng.choice(idx_all, size=n_keep, replace=False) if n_keep > 0 else np.array([], dtype=int)

            Xs = X.flat[chosen]
            Ys = Y.flat[chosen]
            Us = U.flat[chosen]
            Vs = V.flat[chosen]

            # ベクトルの向きはx–yのみ、スケーリング用の大きさも x–y のみで正規化
            mag_xy = np.hypot(Us, Vs)
            nz = mag_xy > 0
            Xs = Xs[nz]
            Ys = Ys[nz]
            Us = (Us[nz] / mag_xy[nz]) * float(vector_size)
            Vs = (Vs[nz] / mag_xy[nz]) * float(vector_size)

            ax.quiver(
                Xs, Ys, Us, Vs,
                color='black', angles='xy', scale_units='xy', scale=1.0, pivot='mid', zorder=4
            )
    elif icon_vector == 2:
        # NaNをゼロに置き換え（streamplotはNaNを処理できないため）
        U_plot = np.nan_to_num(U, nan=0.0)
        V_plot = np.nan_to_num(V, nan=0.0)

        ax.streamplot(
            X, Y, U_plot, V_plot,
            color='black',
            density=2.0,      # 線の密度
            linewidth=0.5,    # 線の太さ
            arrowsize=0.5,    # 矢印サイズ
            zorder=4
        )

    # 軸体裁
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
