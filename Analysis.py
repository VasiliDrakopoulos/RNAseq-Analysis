import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import mygene
import warnings
import os
from scipy.cluster.hierarchy import linkage, leaves_list
warnings.filterwarnings('ignore', category=pd.errors.SettingWithCopyWarning)


# Helper: log2(CPM+1) normalisation
def log2_cpm(counts_df):
    counts = counts_df.copy()
    gene_col = None
    for col in ['Geneid', 'gene_id', 'Gene', 'gene']:
        if col in counts.columns:
            gene_col = col
            break
    if gene_col is None and len(counts.columns) > 0:
        first_col = counts.columns[0]
        if not pd.api.types.is_numeric_dtype(counts[first_col]):
            gene_col = first_col
    if gene_col is not None:
        counts.set_index(gene_col, inplace=True)
    counts = counts.select_dtypes(include='number')
    if counts.empty:
        raise ValueError("No numeric columns found in counts data.")
    lib_sizes = counts.sum(axis=0)
    lib_sizes = lib_sizes.replace(0, np.nan)
    cpm = counts.div(lib_sizes, axis=1) * 1e6
    norm = np.log2(cpm + 1)
    return norm

# Safe conversion helpers
def safe_float(var, default):
    try:
        val = var.get().strip()
        return float(val) if val else float(default)
    except:
        return float(default)

def safe_int(var, default):
    try:
        val = var.get().strip()
        return int(float(val)) if val else int(default)
    except:
        return int(default)

def safe_intvar(var, default):
    try:
        val = var.get()
        return int(val) if val != '' and val is not None else default
    except:
        return default

# Main GUI class
class RNAAnalysisGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("RNA-seq Analysis Figure Generator")
        self.root.geometry("1400x800")
        self.results_df = None
        self.counts_raw = None
        self.counts_norm = None
        self.metadata_df = None
        self.top_genes = None
        self.ensembl_to_symbol = {}
        self.symbol_to_ensembl = {}
        self.ensembl_to_entrez = {}
        self.highlight_file = None
        self.highlight_df = None
        self.group_file = None
        self.group_df = None
        self.create_widgets()

    def create_widgets(self):
        menubar = tk.Menu(self.root)
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Load DESeq2 Results (CSV)", command=self.load_results)
        file_menu.add_command(label="Load Raw Counts (CSV)", command=self.load_counts)
        file_menu.add_command(label="Load Sample Metadata (CSV)", command=self.load_metadata)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)
        menubar.add_cascade(label="File", menu=file_menu)
        self.root.config(menu=menubar)

        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        self.tab_volcano = ttk.Frame(self.notebook)
        self.tab_heatmap = ttk.Frame(self.notebook)
        self.tab_about = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_volcano, text="Volcano Plot")
        self.notebook.add(self.tab_heatmap, text="Heatmap")
        self.notebook.add(self.tab_about, text="About")

        self.build_volcano_tab()
        self.build_heatmap_tab()
        self.build_about_tab()

    # Volcano Tab
    def build_volcano_tab(self):
        main_pane = ttk.PanedWindow(self.tab_volcano, orient=tk.HORIZONTAL)
        main_pane.pack(fill=tk.BOTH, expand=True)
        ctrl_frame = tk.Frame(main_pane)
        main_pane.add(ctrl_frame, weight=1)
        control_canvas = tk.Canvas(ctrl_frame, width=320)
        scrollbar = ttk.Scrollbar(ctrl_frame, orient="vertical", command=control_canvas.yview)
        scrollable_frame = ttk.Frame(control_canvas)
        scrollable_frame.bind("<Configure>", lambda e: control_canvas.configure(scrollregion=control_canvas.bbox("all")))
        control_canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        control_canvas.configure(yscrollcommand=scrollbar.set)
        control_canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self._build_volcano_controls(scrollable_frame)
        self.volcano_canvas_frame = tk.Frame(main_pane)
        main_pane.add(self.volcano_canvas_frame, weight=3)

    def _build_volcano_controls(self, parent):
        ctrl = parent
        tk.Label(ctrl, text="Volcano Plot Options", font=('Arial', 12, 'bold')).pack(pady=5)

        self.enable_secondary = tk.BooleanVar(value=True)
        tk.Checkbutton(ctrl, text="Enable secondary tier", variable=self.enable_secondary).pack(anchor='w', pady=(5,0))

        tk.Label(ctrl, text="Primary P-value cutoff:").pack(anchor='w', pady=(5,0))
        self.pval_cutoff = tk.StringVar(value="0.05")
        tk.Entry(ctrl, textvariable=self.pval_cutoff, width=8).pack(anchor='w')

        tk.Label(ctrl, text="Primary Log2FC cutoff:").pack(anchor='w', pady=(5,0))
        self.lfc_cutoff = tk.StringVar(value="1.0")
        tk.Entry(ctrl, textvariable=self.lfc_cutoff, width=8).pack(anchor='w')

        tk.Label(ctrl, text="Secondary P-value cutoff:").pack(anchor='w', pady=(5,0))
        self.pval_cutoff2 = tk.StringVar(value="0.01")
        tk.Entry(ctrl, textvariable=self.pval_cutoff2, width=8).pack(anchor='w')

        tk.Label(ctrl, text="Secondary Log2FC cutoff:").pack(anchor='w', pady=(5,0))
        self.lfc_cutoff2 = tk.StringVar(value="2.0")
        tk.Entry(ctrl, textvariable=self.lfc_cutoff2, width=8).pack(anchor='w')

        tk.Label(ctrl, text="Colors (NS, Up1, Down1, Up2, Down2):").pack(anchor='w', pady=(5,0))
        self.volcano_colors = tk.StringVar(value="grey,lightcoral,lightblue,red,blue")
        tk.Entry(ctrl, textvariable=self.volcano_colors, width=30).pack(anchor='w')

        tk.Label(ctrl, text="Highlight genes (ID:colour, ID2:colour2):").pack(anchor='w', pady=(5,0))
        self.volcano_highlight_colors = tk.StringVar(value="")
        tk.Entry(ctrl, textvariable=self.volcano_highlight_colors, width=40).pack(anchor='w')

        tk.Label(ctrl, text="Load highlight file (CSV):").pack(anchor='w', pady=(5,0))
        hl_frame = tk.Frame(ctrl)
        hl_frame.pack(anchor='w', fill=tk.X, pady=(0,5))
        self.highlight_file_label = tk.Label(hl_frame, text="No file loaded", relief='sunken', width=30)
        self.highlight_file_label.pack(side=tk.LEFT)
        tk.Button(hl_frame, text="Browse", command=self.load_highlight_file, bg='lightgray').pack(side=tk.LEFT, padx=5)

        tk.Label(ctrl, text="Custom legend labels (colour:label, colour2:label2):").pack(anchor='w', pady=(5,0))
        self.volcano_highlight_legend = tk.StringVar(value="")
        tk.Entry(ctrl, textvariable=self.volcano_highlight_legend, width=40).pack(anchor='w')

        tk.Label(ctrl, text="Legend order (comma-separated colours or labels):").pack(anchor='w', pady=(5,0))
        self.volcano_legend_order = tk.StringVar(value="")
        tk.Entry(ctrl, textvariable=self.volcano_legend_order, width=40).pack(anchor='w')

        tk.Label(ctrl, text="Custom category labels (legend):\n(leave blank to hide)").pack(anchor='w', pady=(5,0))
        tk.Label(ctrl, text="NS:").pack(anchor='w')
        self.volcano_ns_label = tk.StringVar(value="NS")
        tk.Entry(ctrl, textvariable=self.volcano_ns_label, width=15).pack(anchor='w')
        tk.Label(ctrl, text="Up1 (primary up):").pack(anchor='w')
        self.volcano_up1_label = tk.StringVar(value="Up (primary)")
        tk.Entry(ctrl, textvariable=self.volcano_up1_label, width=15).pack(anchor='w')
        tk.Label(ctrl, text="Down1 (primary down):").pack(anchor='w')
        self.volcano_down1_label = tk.StringVar(value="Down (primary)")
        tk.Entry(ctrl, textvariable=self.volcano_down1_label, width=15).pack(anchor='w')
        tk.Label(ctrl, text="Up2 (secondary up):").pack(anchor='w')
        self.volcano_up2_label = tk.StringVar(value="Up (secondary)")
        tk.Entry(ctrl, textvariable=self.volcano_up2_label, width=15).pack(anchor='w')
        tk.Label(ctrl, text="Down2 (secondary down):").pack(anchor='w')
        self.volcano_down2_label = tk.StringVar(value="Down (secondary)")
        tk.Entry(ctrl, textvariable=self.volcano_down2_label, width=15).pack(anchor='w')

        tk.Label(ctrl, text="Rename genes (ID:NewName, ID2:NewName2):").pack(anchor='w', pady=(5,0))
        self.volcano_rename_genes = tk.StringVar(value="")
        tk.Entry(ctrl, textvariable=self.volcano_rename_genes, width=40).pack(anchor='w')

        tk.Label(ctrl, text="Number of top genes to label (0 = none):").pack(anchor='w', pady=(5,0))
        self.volcano_label_n = tk.StringVar(value="10")
        tk.Entry(ctrl, textvariable=self.volcano_label_n, width=6).pack(anchor='w')

        tk.Label(ctrl, text="Repel expansion factor (1.0=default):").pack(anchor='w', pady=(5,0))
        self.repel_expand = tk.StringVar(value="1.0")
        tk.Entry(ctrl, textvariable=self.repel_expand, width=6).pack(anchor='w')

        tk.Label(ctrl, text="Label offset (x, y) after repel (points):").pack(anchor='w', pady=(5,0))
        frame_offset = tk.Frame(ctrl)
        frame_offset.pack(anchor='w', pady=(0,5))
        self.volcano_label_offset_x = tk.StringVar(value="0")
        self.volcano_label_offset_y = tk.StringVar(value="0")
        tk.Entry(frame_offset, textvariable=self.volcano_label_offset_x, width=4).pack(side=tk.LEFT)
        tk.Label(frame_offset, text=",").pack(side=tk.LEFT)
        tk.Entry(frame_offset, textvariable=self.volcano_label_offset_y, width=4).pack(side=tk.LEFT)

        tk.Label(ctrl, text="Title:").pack(anchor='w', pady=(5,0))
        self.volcano_title = tk.StringVar(value="Volcano Plot")
        tk.Entry(ctrl, textvariable=self.volcano_title, width=25).pack(anchor='w')

        tk.Label(ctrl, text="Title font size:").pack(anchor='w', pady=(5,0))
        self.volcano_title_font = tk.StringVar(value="14")
        tk.Entry(ctrl, textvariable=self.volcano_title_font, width=6).pack(anchor='w')

        tk.Label(ctrl, text="Axis label font size:").pack(anchor='w', pady=(5,0))
        self.volcano_axis_label_font = tk.StringVar(value="12")
        tk.Entry(ctrl, textvariable=self.volcano_axis_label_font, width=6).pack(anchor='w')

        tk.Label(ctrl, text="Tick label font size:").pack(anchor='w', pady=(5,0))
        self.volcano_tick_font = tk.StringVar(value="10")
        tk.Entry(ctrl, textvariable=self.volcano_tick_font, width=6).pack(anchor='w')

        tk.Label(ctrl, text="Figure width, height:").pack(anchor='w', pady=(5,0))
        self.volcano_figsize = tk.StringVar(value="10,7")
        tk.Entry(ctrl, textvariable=self.volcano_figsize, width=10).pack(anchor='w')

        tk.Button(ctrl, text="Generate Volcano", command=self.plot_volcano, bg='lightblue', font=('Arial', 10, 'bold')).pack(pady=10)
        tk.Button(ctrl, text="Save Figure", command=lambda: self.save_figure("volcano"), bg='lightgreen', font=('Arial', 10, 'bold')).pack(pady=5)

    def _resolve_gene_id(self, identifier):
        if identifier.startswith('ENSMUSG'):
            return identifier
        if self.symbol_to_ensembl and identifier in self.symbol_to_ensembl:
            return self.symbol_to_ensembl[identifier]
        return identifier

    def load_highlight_file(self):
        fname = filedialog.askopenfilename(parent=self.root, title="Select Highlight CSV", filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if fname:
            try:
                df = pd.read_csv(fname)
                if 'Gene' not in df.columns or 'Colour' not in df.columns:
                    messagebox.showerror("Error", "CSV must have 'Gene' and 'Colour' columns.")
                    return
                self.highlight_file = fname
                self.highlight_df = df
                self.highlight_file_label.config(text=os.path.basename(fname))
                messagebox.showinfo("Success", f"Loaded {len(df)} highlight entries.")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file:\n{str(e)}")

    def _get_highlight_dict_from_file(self):
        if self.highlight_df is None:
            return {}, {}
        hl_dict, leg_dict = {}, {}
        for _, row in self.highlight_df.iterrows():
            gene = str(row['Gene']).strip()
            colour = str(row['Colour']).strip()
            if gene and colour:
                gene_id = self._resolve_gene_id(gene)
                hl_dict[gene_id] = colour
                if 'Label' in row and pd.notna(row['Label']):
                    label = str(row['Label']).strip()
                    if label:
                        leg_dict[colour] = label
        return hl_dict, leg_dict

    def plot_volcano(self):
        try:
            if self.results_df is None:
                messagebox.showerror("Error", "Load DESeq2 results first.")
                return
            df = self.results_df.copy()
            df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
            df['log2FoldChange'] = pd.to_numeric(df['log2FoldChange'], errors='coerce')
            df = df.dropna(subset=['padj', 'log2FoldChange'])
            if df.empty:
                messagebox.showerror("Error", "No valid numeric data after conversion.")
                return

            p_cut1 = safe_float(self.pval_cutoff, 0.05)
            lfc_cut1 = safe_float(self.lfc_cutoff, 1.0)
            p_cut2 = safe_float(self.pval_cutoff2, 0.01)
            lfc_cut2 = safe_float(self.lfc_cutoff2, 2.0)
            title = self.volcano_title.get().strip()
            use_secondary = self.enable_secondary.get()
            expand_factor = safe_float(self.repel_expand, 1.0)
            offset_x = float(safe_float(self.volcano_label_offset_x, 0))
            offset_y = float(safe_float(self.volcano_label_offset_y, 0))

            colors_raw = self.volcano_colors.get().strip()
            if colors_raw:
                parts = [p.strip() for p in colors_raw.replace(',', ' ').split() if p.strip()]
            else:
                parts = []
            default_colors = ['grey', 'lightcoral', 'lightblue', 'red', 'blue']
            colors = (parts + default_colors)[:5]

            title_fs = safe_int(self.volcano_title_font, 14)
            label_fs = safe_int(self.volcano_axis_label_font, 12)
            tick_fs = safe_int(self.volcano_tick_font, 10)
            try:
                w, h = map(float, self.volcano_figsize.get().split(','))
            except:
                w, h = 10, 7

            ns_label = self.volcano_ns_label.get().strip()
            up1_label = self.volcano_up1_label.get().strip()
            down1_label = self.volcano_down1_label.get().strip()
            up2_label = self.volcano_up2_label.get().strip()
            down2_label = self.volcano_down2_label.get().strip()

            rename_dict = {}
            rename_str = self.volcano_rename_genes.get().strip()
            if rename_str:
                for item in rename_str.split(','):
                    if ':' in item:
                        k, v = item.split(':', 1)
                        rename_dict[self._resolve_gene_id(k.strip())] = v.strip()

            highlight_dict, highlight_legend = self._get_highlight_dict_from_file()
            hl_str = self.volcano_highlight_colors.get().strip()
            if hl_str:
                for item in hl_str.split(','):
                    if ':' in item:
                        k, v = item.split(':', 1)
                        highlight_dict[self._resolve_gene_id(k.strip())] = v.strip()
            leg_str = self.volcano_highlight_legend.get().strip()
            if leg_str:
                for item in leg_str.split(','):
                    if ':' in item:
                        k, v = item.split(':', 1)
                        highlight_legend[k.strip()] = v.strip()

            if not use_secondary:
                ns = (df['padj'] >= p_cut1) | (abs(df['log2FoldChange']) < lfc_cut1)
                up = (df['padj'] < p_cut1) & (df['log2FoldChange'] > lfc_cut1)
                down = (df['padj'] < p_cut1) & (df['log2FoldChange'] < -lfc_cut1)
                df['sig'] = 'NS'
                df.loc[up, 'sig'] = 'Up'
                df.loc[down, 'sig'] = 'Down'
                display_map = {'NS': ns_label, 'Up': up1_label, 'Down': down1_label}
                color_map = {'NS': colors[0], 'Up': colors[1], 'Down': colors[2]}
                categories = ['NS', 'Up', 'Down']
            else:
                ns = (df['padj'] >= p_cut1) | (abs(df['log2FoldChange']) < lfc_cut1)
                up2 = (df['padj'] < p_cut2) & (df['log2FoldChange'] > lfc_cut2)
                down2 = (df['padj'] < p_cut2) & (df['log2FoldChange'] < -lfc_cut2)
                up1 = (df['padj'] < p_cut1) & (df['log2FoldChange'] > lfc_cut1) & ~up2
                down1 = (df['padj'] < p_cut1) & (df['log2FoldChange'] < -lfc_cut1) & ~down2
                df['sig'] = 'NS'
                df.loc[up1, 'sig'] = 'Up1'
                df.loc[down1, 'sig'] = 'Down1'
                df.loc[up2, 'sig'] = 'Up2'
                df.loc[down2, 'sig'] = 'Down2'
                display_map = {'NS': ns_label, 'Up1': up1_label, 'Down1': down1_label, 'Up2': up2_label, 'Down2': down2_label}
                color_map = {'NS': colors[0], 'Up1': colors[1], 'Down1': colors[2], 'Up2': colors[3], 'Down2': colors[4]}
                categories = ['NS', 'Up1', 'Down1', 'Up2', 'Down2']

            fig, ax = plt.subplots(figsize=(w, h))
            hl_points = {}
            for name in categories:
                label_str = display_map[name]
                if label_str:
                    if name == 'NS':
                        full_label = label_str
                    elif name.startswith('Up'):
                        full_label = f"{label_str} (p<{p_cut2 if '2' in name else p_cut1}, FC>{lfc_cut2 if '2' in name else lfc_cut1})"
                    elif name.startswith('Down'):
                        full_label = f"{label_str} (p<{p_cut2 if '2' in name else p_cut1}, FC<{-lfc_cut2 if '2' in name else -lfc_cut1})"
                    else:
                        full_label = label_str
                else:
                    full_label = None

                group = df[df['sig'] == name]
                if group.empty:
                    if full_label is not None:
                        ax.scatter([], [], c=color_map[name], label=full_label, s=0, alpha=0.6)
                    continue

                base_color = color_map[name]
                non_hl = []
                for _, row in group.iterrows():
                    gene_id = row['Geneid']
                    if gene_id in highlight_dict:
                        hl_color = highlight_dict[gene_id]
                        hl_points.setdefault(hl_color, {'x': [], 'y': []})
                        hl_points[hl_color]['x'].append(row['log2FoldChange'])
                        hl_points[hl_color]['y'].append(-np.log10(row['padj']))
                    else:
                        non_hl.append(row)
                if non_hl:
                    non_hl_df = pd.DataFrame(non_hl)
                    ax.scatter(non_hl_df['log2FoldChange'], -np.log10(non_hl_df['padj']),
                               c=base_color, label=full_label, alpha=0.6,
                               s=20 if name=='NS' else (40 if name in ['Up2','Down2'] else 30))
                else:
                    if full_label is not None:
                        ax.scatter([], [], c=base_color, label=full_label, s=0, alpha=0.6)

            order_str = self.volcano_legend_order.get().strip()
            if order_str:
                user_order = [x.strip() for x in order_str.split(',') if x.strip()]
                def sort_key(item):
                    colour = item[0]
                    label = highlight_legend.get(colour, colour)
                    if label in user_order:
                        return user_order.index(label)
                    if colour in user_order:
                        return user_order.index(colour)
                    return len(user_order)
            else:
                def sort_key(item):
                    return highlight_legend.get(item[0], item[0])

            for hl_color, points in sorted(hl_points.items(), key=sort_key):
                ax.scatter(points['x'], points['y'], c=hl_color,
                           label=highlight_legend.get(hl_color, None),
                           alpha=0.8, edgecolors='black', linewidth=0.5, s=50)

            if use_secondary:
                ax.axhline(-np.log10(p_cut1), linestyle='--', color='k', alpha=0.3)
                ax.axhline(-np.log10(p_cut2), linestyle=':', color='k', alpha=0.3)
                ax.axvline(-lfc_cut1, linestyle='--', color='k', alpha=0.3)
                ax.axvline(lfc_cut1, linestyle='--', color='k', alpha=0.3)
                ax.axvline(-lfc_cut2, linestyle=':', color='k', alpha=0.3)
                ax.axvline(lfc_cut2, linestyle=':', color='k', alpha=0.3)
            else:
                ax.axhline(-np.log10(p_cut1), linestyle='--', color='k', alpha=0.3)
                ax.axvline(-lfc_cut1, linestyle='--', color='k', alpha=0.3)
                ax.axvline(lfc_cut1, linestyle='--', color='k', alpha=0.3)

            ax.set_xlabel('log2 Fold Change', fontsize=label_fs)
            ax.set_ylabel('-log10 adjusted p-value', fontsize=label_fs)
            if title:
                ax.set_title(title, fontsize=title_fs)
            ax.tick_params(labelsize=tick_fs)
            ax.legend(fontsize=tick_fs-1, loc='best')
            sns.despine()

            label_str = self.volcano_label_n.get().strip()
            label_n = int(float(label_str)) if label_str else 0
            if label_n > 0 and self.ensembl_to_symbol:
                top_genes = df.nsmallest(label_n, 'padj')
                try:
                    from adjustText import adjust_text
                    texts = []
                    for _, row in top_genes.iterrows():
                        gene_id = row['Geneid']
                        symbol = rename_dict.get(gene_id, self.ensembl_to_symbol.get(gene_id, gene_id))
                        if symbol.startswith('ENSMUSG'):
                            symbol = self.ensembl_to_symbol.get(gene_id, gene_id)
                        txt = ax.text(row['log2FoldChange'], -np.log10(row['padj']), symbol,
                                      fontsize=tick_fs-1, alpha=0.8,
                                      bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3))
                        texts.append(txt)
                    adjust_text(texts, ax=ax, expand_points=(expand_factor, expand_factor),
                                arrowprops=dict(arrowstyle='->', color='red', lw=0.5))
                    for txt in texts:
                        x, y = txt.get_position()
                        txt.set_position((x + offset_x, y + offset_y))
                except ImportError:
                    for _, row in top_genes.iterrows():
                        gene_id = row['Geneid']
                        symbol = rename_dict.get(gene_id, self.ensembl_to_symbol.get(gene_id, gene_id))
                        if symbol.startswith('ENSMUSG'):
                            symbol = self.ensembl_to_symbol.get(gene_id, gene_id)
                        ax.annotate(symbol, (row['log2FoldChange'], -np.log10(row['padj'])),
                                    fontsize=tick_fs-1, alpha=0.8,
                                    xytext=(offset_x, offset_y), textcoords='offset points',
                                    bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3))

            plt.tight_layout()
            self.display_figure(fig, self.volcano_canvas_frame)
        except Exception as e:
            messagebox.showerror("Volcano Error", f"An error occurred:\n{str(e)}")
            import traceback
            traceback.print_exc()

    # Heatmap Tab
    def build_heatmap_tab(self):
        main_pane = ttk.PanedWindow(self.tab_heatmap, orient=tk.HORIZONTAL)
        main_pane.pack(fill=tk.BOTH, expand=True)
        ctrl = tk.Frame(main_pane)
        main_pane.add(ctrl, weight=1)

        tk.Label(ctrl, text="Heatmap Options", font=('Arial', 12, 'bold')).pack(pady=5)

        tk.Label(ctrl, text="Number of top genes (single mode):").pack(anchor='w')
        self.top_n = tk.IntVar(value=10)
        tk.Entry(ctrl, textvariable=self.top_n, width=8).pack(anchor='w')

        tk.Label(ctrl, text="Sort by (single mode):").pack(anchor='w')
        self.sort_by = tk.StringVar(value="padj")
        ttk.Combobox(ctrl, textvariable=self.sort_by, values=["padj", "log2FoldChange"], state='readonly', width=15).pack(anchor='w')

        tk.Label(ctrl, text="Color palette (single mode):").pack(anchor='w')
        self.heatmap_cmap = tk.StringVar(value="RdBu_r")
        ttk.Combobox(ctrl, textvariable=self.heatmap_cmap, values=["RdBu_r", "viridis", "plasma", "coolwarm", "seismic"], state='readonly', width=15).pack(anchor='w')

        tk.Label(ctrl, text="Colormaps (comma-separated, for multi-block):").pack(anchor='w')
        self.heatmap_cmaps = tk.StringVar(value="viridis,plasma,inferno,magma,cividis")
        tk.Entry(ctrl, textvariable=self.heatmap_cmaps, width=40).pack(anchor='w')

        tk.Label(ctrl, text="Title:").pack(anchor='w')
        self.heatmap_title = tk.StringVar(value="Top Genes Heatmap")
        tk.Entry(ctrl, textvariable=self.heatmap_title, width=20).pack(anchor='w')

        tk.Label(ctrl, text="Rename genes (ID:NewName, ID2:NewName2):").pack(anchor='w', pady=(5,0))
        self.heatmap_rename_genes = tk.StringVar(value="")
        tk.Entry(ctrl, textvariable=self.heatmap_rename_genes, width=40).pack(anchor='w')

        tk.Label(ctrl, text="Custom sample labels (comma-separated):").pack(anchor='w')
        self.heatmap_labels = tk.StringVar(value="")
        tk.Entry(ctrl, textvariable=self.heatmap_labels, width=30).pack(anchor='w')

        tk.Label(ctrl, text="Metadata sample column (must match counts headers):").pack(anchor='w', pady=(5,0))
        self.sample_col_var = tk.StringVar()
        self.sample_col_dropdown = ttk.Combobox(ctrl, textvariable=self.sample_col_var, state='readonly', width=20)
        self.sample_col_dropdown.pack(anchor='w', pady=(0,5))

        tk.Label(ctrl, text="Group annotation column (e.g., 'condition'):").pack(anchor='w')
        self.group_col_var = tk.StringVar()
        self.group_col_dropdown = ttk.Combobox(ctrl, textvariable=self.group_col_var, state='readonly', width=20)
        self.group_col_dropdown.pack(anchor='w', pady=(0,5))

        tk.Button(ctrl, text="Refresh metadata columns", command=self.update_metadata_dropdowns, bg='lightgray').pack(anchor='w', pady=5)

        self.heatmap_cluster_cols = tk.BooleanVar(value=True)
        tk.Checkbutton(ctrl, text="Cluster columns?", variable=self.heatmap_cluster_cols).pack(anchor='w')

        tk.Label(ctrl, text="Title font size:").pack(anchor='w')
        self.heatmap_title_font = tk.IntVar(value=14)
        tk.Entry(ctrl, textvariable=self.heatmap_title_font, width=6).pack(anchor='w')

        tk.Label(ctrl, text="Tick label font size:").pack(anchor='w')
        self.heatmap_tick_font = tk.IntVar(value=10)
        tk.Entry(ctrl, textvariable=self.heatmap_tick_font, width=6).pack(anchor='w')

        tk.Label(ctrl, text="Figure width, height:").pack(anchor='w')
        self.heatmap_figsize = tk.StringVar(value="10,6")
        tk.Entry(ctrl, textvariable=self.heatmap_figsize, width=10).pack(anchor='w')

        tk.Label(ctrl, text="Load group file (CSV) for multi‑block heatmap:").pack(anchor='w', pady=(5,0))
        group_frame = tk.Frame(ctrl)
        group_frame.pack(anchor='w', fill=tk.X, pady=(0,5))
        self.group_file_label = tk.Label(group_frame, text="No file loaded", relief='sunken', width=30)
        self.group_file_label.pack(side=tk.LEFT)
        tk.Button(group_frame, text="Browse", command=self.load_group_file, bg='lightgray').pack(side=tk.LEFT, padx=5)
        tk.Label(ctrl, text="(If loaded, single‑mode options above are ignored)", font=('Arial', 8)).pack(anchor='w')

        tk.Button(ctrl, text="Generate Heatmap", command=self.plot_heatmap, bg='lightblue').pack(pady=10)
        tk.Button(ctrl, text="Save Figure", command=lambda: self.save_figure("heatmap"), bg='lightgreen').pack(pady=5)
        self.heatmap_canvas_frame = tk.Frame(main_pane)
        main_pane.add(self.heatmap_canvas_frame, weight=3)

    def update_metadata_dropdowns(self):
        if self.metadata_df is None:
            messagebox.showwarning("Warning", "No metadata loaded. Load a metadata file first.")
            return
        cols = list(self.metadata_df.columns)
        self.sample_col_dropdown['values'] = cols
        self.group_col_dropdown['values'] = cols
        for col in cols:
            if col.lower() in ['sample', 'sampleid', 'id']:
                self.sample_col_var.set(col)
            if col.lower() in ['condition', 'group', 'treatment', 'sex']:
                self.group_col_var.set(col)
        messagebox.showinfo("Info", f"Metadata columns: {', '.join(cols)}")

    def load_group_file(self):
        fname = filedialog.askopenfilename(parent=self.root, title="Select Group CSV", filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if fname:
            try:
                df = pd.read_csv(fname)
                if 'Gene' not in df.columns or 'Group' not in df.columns:
                    messagebox.showerror("Error", "CSV must have 'Gene' and 'Group' columns.")
                    return
                self.group_file = fname
                self.group_df = df
                self.group_file_label.config(text=os.path.basename(fname))
                messagebox.showinfo("Success", f"Loaded {len(df)} genes in groups.")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file:\n{str(e)}")

    def _cluster_columns(self, mat_scaled, group_info=None):
        if group_info is None:
            Z = linkage(mat_scaled.T, method='average', metric='euclidean')
            order = leaves_list(Z)
            return mat_scaled.iloc[:, order], mat_scaled.columns[order].tolist()
        else:
            col_names = mat_scaled.columns.tolist()
            new_order = []
            for start, end, _ in group_info:
                subset_cols = col_names[start:end+1]
                if len(subset_cols) <= 1:
                    new_order.extend(subset_cols)
                else:
                    submat = mat_scaled[subset_cols]
                    Z = linkage(submat.T, method='average', metric='euclidean')
                    leaves = leaves_list(Z)
                    new_order.extend([subset_cols[i] for i in leaves])
            mat_reordered = mat_scaled[new_order]
            return mat_reordered, new_order

    def _build_group_info(self, sample_col, group_col, col_order):
        if not sample_col or not group_col or self.metadata_df is None:
            return None
        meta = self.metadata_df.copy()
        if sample_col not in meta.columns or group_col not in meta.columns:
            return None
        meta = meta.set_index(sample_col)
        meta.index = meta.index.astype(str).str.strip()
        col_order = [c.strip() for c in col_order]
        common = [s for s in col_order if s in meta.index]
        if not common:
            return None
        groups = meta.loc[common, group_col]
        group_info = []
        current = groups.iloc[0]
        start = 0
        for i, g in enumerate(groups):
            if g != current:
                group_info.append((start, i-1, current))
                start = i
                current = g
        group_info.append((start, len(groups)-1, current))
        return group_info

    def plot_heatmap(self):
        try:
            if self.results_df is None:
                messagebox.showerror("Error", "Load DESeq2 results first.")
                return
            if self.counts_raw is None:
                messagebox.showerror("Error", "Load raw counts file first.")
                return
            if self.counts_norm is None:
                self.counts_norm = log2_cpm(self.counts_raw)

            # Multi‑block heatmap (group file)
            if self.group_df is not None:
                group_map = {}
                groups_order = []
                for _, row in self.group_df.iterrows():
                    gene = str(row['Gene']).strip()
                    group = str(row['Group']).strip()
                    if gene and group:
                        gene_id = self._resolve_gene_id(gene)
                        group_map[gene_id] = group
                        if group not in groups_order:
                            groups_order.append(group)

                df = self.results_df.copy()
                df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
                df['log2FoldChange'] = pd.to_numeric(df['log2FoldChange'], errors='coerce')
                df = df.dropna(subset=['padj', 'log2FoldChange'])
                df['group'] = df['Geneid'].map(group_map)
                df = df.dropna(subset=['group'])
                if df.empty:
                    messagebox.showerror("Error", "None of the genes in the group file are present in the DESeq2 results.")
                    return

                norm = self.counts_norm
                present_genes = [g for g in df['Geneid'].tolist() if g in norm.index]
                if not present_genes:
                    messagebox.showerror("Error", "None of the genes found in counts. Check Geneid matching.")
                    return
                mat = norm.loc[present_genes]

                group_data = {}
                sort_col = self.sort_by.get()
                for group in groups_order:
                    genes_in_group = df[df['group'] == group]['Geneid'].tolist()
                    if sort_col == 'padj':
                        genes_sorted = df[df['group'] == group].sort_values('padj')['Geneid'].tolist()
                    else:
                        genes_sorted = df[df['group'] == group].sort_values('log2FoldChange', ascending=False)['Geneid'].tolist()
                    mat_sub = mat.loc[[g for g in genes_sorted if g in mat.index]]
                    if mat_sub.empty:
                        continue
                    mat_scaled = (mat_sub.T - mat_sub.mean(axis=1)).T
                    group_data[group] = mat_scaled

                if not group_data:
                    messagebox.showerror("Error", "No valid groups with genes found.")
                    return

                custom_labels_str = self.heatmap_labels.get().strip()
                if custom_labels_str:
                    labels = [x.strip() for x in custom_labels_str.split(',')]
                    if len(labels) == mat.shape[1]:
                        col_labels = labels
                    else:
                        messagebox.showwarning("Warning", f"Custom label count mismatch. Using sample names.")
                        col_labels = mat.columns.tolist()
                else:
                    col_labels = mat.columns.tolist()

                try:
                    w, h = map(float, self.heatmap_figsize.get().split(','))
                except:
                    w, h = 10, 6
                total_genes = sum([len(g) for g in group_data.values()])
                h = max(h, float(total_genes * 0.3 + 2))

                n_groups = len(group_data)
                fig, axes = plt.subplots(n_groups, 1, figsize=(w, h), sharex=True)
                if n_groups == 1:
                    axes = [axes]

                cmaps_str = self.heatmap_cmaps.get().strip()
                if cmaps_str:
                    cmaps = [c.strip() for c in cmaps_str.split(',') if c.strip()]
                else:
                    cmaps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis']
                if len(cmaps) < n_groups:
                    cmaps = (cmaps * (n_groups // len(cmaps) + 1))[:n_groups]
                else:
                    cmaps = cmaps[:n_groups]

                tick_fs = safe_intvar(self.heatmap_tick_font, 10)
                rename_dict = {}
                rename_str = self.heatmap_rename_genes.get().strip()
                if rename_str:
                    for item in rename_str.split(','):
                        if ':' in item:
                            k, v = item.split(':', 1)
                            rename_dict[self._resolve_gene_id(k.strip())] = v.strip()

                sample_col = self.sample_col_var.get().strip()
                group_col = self.group_col_var.get().strip()
                group_info = self._build_group_info(sample_col, group_col, mat.columns.tolist())
                n_cols = mat.shape[1]

                for ax, (group_name, mat_scaled) in zip(axes, group_data.items()):
                    y_labels = []
                    for g in mat_scaled.index:
                        symbol = rename_dict.get(g, self.ensembl_to_symbol.get(g, g))
                        if symbol.startswith('ENSMUSG'):
                            symbol = self.ensembl_to_symbol.get(g, g)
                        y_labels.append(symbol)
                    cmap_idx = list(group_data.keys()).index(group_name)
                    cmap = cmaps[cmap_idx]

                    sns.heatmap(mat_scaled, cmap=cmap, cbar=True, ax=ax,
                                xticklabels=False,
                                yticklabels=y_labels,
                                cbar_kws={'label': f'Z-score ({group_name})'})
                    ax.set_ylabel(group_name, fontsize=tick_fs+2, weight='bold')
                    ax.tick_params(labelsize=tick_fs)
                    if ax != axes[-1]:
                        ax.set_xlabel('')
                    ax.set_xlim(-0.5, n_cols - 0.5)

                bottom_ax = axes[-1]
                if group_info is not None:
                    bottom_ax.set_xticks([])
                    for start, end, label in group_info:
                        mid = (start + end) / 2
                        x_frac = (mid + 0.5) / n_cols
                        bottom_ax.text(x_frac, -0.02, label, ha='center', va='top',
                                       fontsize=tick_fs+1, weight='bold', transform=bottom_ax.transAxes)

                    # Draw dotted lines at group boundaries: offset +1.0 (right edge of last sample)
                    for start, end, _ in group_info[:-1]:
                        x_pos = end + 1.0
                        for ax in axes:
                            ax.axvline(x_pos, color='black', linestyle=':', linewidth=1.5, alpha=0.8, zorder=10)
                    plt.subplots_adjust(bottom=0.08)
                else:
                    bottom_ax.set_xticks(np.arange(n_cols) + 0.5)
                    bottom_ax.set_xticklabels(col_labels, rotation=45, ha='right', fontsize=tick_fs)

                fig.suptitle(self.heatmap_title.get(), fontsize=safe_intvar(self.heatmap_title_font, 14))
                plt.tight_layout()
                self.display_figure(fig, self.heatmap_canvas_frame)
                return

            # Single‑block heatmap
            df = self.results_df.dropna(subset=['padj']).copy()
            df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
            df = df.dropna(subset=['padj'])
            n = safe_intvar(self.top_n, 10)
            sort_col = self.sort_by.get()
            if sort_col == 'padj':
                df_sorted = df.sort_values('padj').head(n)
            else:
                df_sorted = df.sort_values('log2FoldChange', ascending=False).head(n)

            gene_list = df_sorted['Geneid'].tolist()
            self.top_genes = gene_list

            norm = self.counts_norm
            present_genes = [g for g in gene_list if g in norm.index]
            if not present_genes:
                messagebox.showerror("Error", "None of the top genes found in counts.")
                return
            mat = norm.loc[present_genes]
            mat_scaled = (mat.T - mat.mean(axis=1)).T

            rename_dict = {}
            rename_str = self.heatmap_rename_genes.get().strip()
            if rename_str:
                for item in rename_str.split(','):
                    if ':' in item:
                        k, v = item.split(':', 1)
                        rename_dict[self._resolve_gene_id(k.strip())] = v.strip()
            y_labels = []
            for g in mat.index:
                symbol = rename_dict.get(g, self.ensembl_to_symbol.get(g, g))
                if symbol.startswith('ENSMUSG'):
                    symbol = self.ensembl_to_symbol.get(g, g)
                y_labels.append(symbol)

            custom_labels_str = self.heatmap_labels.get().strip()
            if custom_labels_str:
                labels = [x.strip() for x in custom_labels_str.split(',')]
                if len(labels) == mat_scaled.shape[1]:
                    col_labels = labels
                else:
                    messagebox.showwarning("Warning", f"Custom label count mismatch. Using sample names.")
                    col_labels = mat_scaled.columns.tolist()
            else:
                col_labels = mat_scaled.columns.tolist()

            try:
                w, h = map(float, self.heatmap_figsize.get().split(','))
            except:
                w, h = 10, 6
            h = max(h, float(n * 0.4 + 2))

            cluster_on = self.heatmap_cluster_cols.get()
            sample_col = self.sample_col_var.get().strip()
            group_col = self.group_col_var.get().strip()
            group_info = self._build_group_info(sample_col, group_col, mat_scaled.columns.tolist())

            if cluster_on and group_info is not None:
                mat_scaled, new_order = self._cluster_columns(mat_scaled, group_info)
                group_info = self._build_group_info(sample_col, group_col, mat_scaled.columns.tolist())
                if custom_labels_str and len(labels) == mat_scaled.shape[1]:
                    col_labels = mat_scaled.columns.tolist()
                else:
                    col_labels = mat_scaled.columns.tolist()

            fig, ax = plt.subplots(figsize=(w, h))
            tick_fs = safe_intvar(self.heatmap_tick_font, 10)

            if group_info is not None:
                sns.heatmap(mat_scaled, cmap=self.heatmap_cmap.get(), cbar=True, ax=ax,
                            xticklabels=False, yticklabels=y_labels,
                            cbar_kws={'label': 'Z-score'})
                ax.set_xticks([])
                n_cols = mat_scaled.shape[1]
                ax.set_xlim(-0.5, n_cols - 0.5)
                for start, end, label in group_info:
                    mid = (start + end) / 2
                    x_frac = (mid + 0.5) / n_cols
                    ax.text(x_frac, -0.02, label, ha='center', va='top',
                            fontsize=tick_fs+1, weight='bold', transform=ax.transAxes)
                for start, end, _ in group_info[:-1]:
                    x_pos = end + 1.0 
                    ax.axvline(x_pos, color='black', linestyle=':', linewidth=1.5, alpha=0.8, zorder=10)
                plt.subplots_adjust(bottom=0.08)
            else:
                sns.heatmap(mat_scaled, cmap=self.heatmap_cmap.get(), cbar=True, ax=ax,
                            xticklabels=col_labels, yticklabels=y_labels,
                            cbar_kws={'label': 'Z-score'})
                ax.set_xticklabels(col_labels, rotation=45, ha='right', fontsize=tick_fs)

            ax.set_title(self.heatmap_title.get(), fontsize=safe_intvar(self.heatmap_title_font, 14))
            ax.tick_params(labelsize=tick_fs)
            ax.set_ylabel('Genes', fontsize=tick_fs+2)

            plt.tight_layout()
            self.display_figure(fig, self.heatmap_canvas_frame)

        except Exception as e:
            messagebox.showerror("Heatmap Error", f"An error occurred:\n{str(e)}")
            import traceback
            traceback.print_exc()

    # About Tab
    def build_about_tab(self):
        frame = self.tab_about
        tk.Label(frame, text="RNA-seq Analysis Figure Generator", font=('Arial', 14, 'bold')).pack(pady=10)
        tk.Label(frame, text="Volcano & Heatmap only", font=('Arial', 10)).pack()
        tk.Label(frame, text="\nFeatures:\n"
                              "- Multi-tier volcano with per-gene colour highlighting\n"
                              "- Heatmap: single or multi‑block (load group CSV)\n"
                              "- Custom colormaps per block\n"
                              "- ID mapping via mygene (automatic)",
                 justify=tk.LEFT, padx=10).pack(anchor='w')

    # Load methods
    def load_results(self):
        fname = filedialog.askopenfilename(parent=self.root, title="Select DESeq2 Results CSV", filetypes=[("CSV files", "*.csv")])
        if fname:
            self.results_df = pd.read_csv(fname)
            self._map_gene_ids_mygene()
            if self.ensembl_to_symbol:
                self.symbol_to_ensembl = {v: k for k, v in self.ensembl_to_symbol.items() if v != k}
                messagebox.showinfo("Success", f"Loaded {len(self.results_df)} genes, mapped {len(self.ensembl_to_symbol)} to symbols.")
            else:
                messagebox.showwarning("Warning", "Gene mapping failed. No symbol mapping available for volcano labels.")

    def _map_gene_ids_mygene(self):
        if self.results_df is None:
            return
        ids = self.results_df['Geneid'].tolist()
        clean = list(set([i.split('.')[0] for i in ids if isinstance(i, str) and i.startswith('ENSMUSG')]))
        if not clean:
            return
        try:
            mg = mygene.MyGeneInfo()
            results = mg.querymany(clean, scopes='ensembl.gene', fields='symbol,entrezgene',
                                   species='mouse', returnall=True, verbose=False, chunk_size=100, sleep=0.5)
            self.ensembl_to_symbol = {}
            self.ensembl_to_entrez = {}
            for item in results['out']:
                if 'query' in item:
                    q = item['query']
                    if 'symbol' in item:
                        self.ensembl_to_symbol[q] = item['symbol']
                    if 'entrezgene' in item:
                        self.ensembl_to_entrez[q] = str(item['entrezgene'])
            for eid in clean:
                self.ensembl_to_symbol.setdefault(eid, eid)
                self.ensembl_to_entrez.setdefault(eid, '')
        except Exception as e:
            print(f"mygene error: {e}")
            self.ensembl_to_symbol = {}

    def load_counts(self):
        fname = filedialog.askopenfilename(parent=self.root, title="Select Raw Counts CSV", filetypes=[("CSV files", "*.csv")])
        if fname:
            df = pd.read_csv(fname)
            if len(df.columns) > 0 and not pd.api.types.is_numeric_dtype(df.iloc[:, 0]):
                df.set_index(df.columns[0], inplace=True)
            self.counts_raw = df
            self.counts_norm = None
            messagebox.showinfo("Success", f"Loaded counts with {len(self.counts_raw.columns)} samples.")

    def load_metadata(self):
        fname = filedialog.askopenfilename(parent=self.root, title="Select Sample Metadata CSV", filetypes=[("CSV files", "*.csv")])
        if fname:
            self.metadata_df = pd.read_csv(fname)
            self.update_metadata_dropdowns()
            messagebox.showinfo("Success", f"Loaded metadata with {len(self.metadata_df)} rows.")

    def display_figure(self, fig, parent_frame):
        for widget in parent_frame.winfo_children():
            widget.destroy()
        canvas = FigureCanvasTkAgg(fig, master=parent_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.current_fig = fig
        self.current_canvas = canvas

    def save_figure(self, plot_type):
        if not hasattr(self, 'current_fig'):
            messagebox.showerror("Error", "Generate a plot first.")
            return
        fname = filedialog.asksaveasfilename(parent=self.root, defaultextension=".png",
                                             filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
        if fname:
            self.current_fig.savefig(fname, dpi=300, bbox_inches='tight')
            messagebox.showinfo("Saved", f"Figure saved to {fname}")

if __name__ == "__main__":
    root = tk.Tk()
    app = RNAAnalysisGUI(root)
    root.mainloop()